function [Hest,posteriorVariance,dataMask] = assemble_burst_channel_estimates_v15_1( ...
    blocks,blockDataIndices,nFFT,nRx,nTx,nSymbols)
%ASSEMBLE_BURST_CHANNEL_ESTIMATES_V15_1 Assemble only per-burst causal CSI.
% Each block's estimate was generated at that block's pilot occasion. Copying
% these estimates by their own data-symbol indices prevents future pilot
% observations from leaking into detection of earlier data bursts.

if ~isstruct(blocks) || isempty(blocks) || ~iscell(blockDataIndices) || ...
        numel(blocks) ~= numel(blockDataIndices)
    error('assemble_burst_channel_estimates_v15_1:InvalidBlocks', ...
        'blocks and blockDataIndices must have the same nonzero number of entries.');
end
counts=[nFFT nRx nTx nSymbols];
if any(~isfinite(counts)) || any(counts<1) || any(counts~=floor(counts))
    error('assemble_burst_channel_estimates_v15_1:InvalidDimensions', ...
        'All channel dimensions must be positive integers.');
end
Hest=complex(zeros(nFFT,nRx,nTx,nSymbols));
posteriorVariance=zeros(nTx,nSymbols);
dataMask=false(1,nSymbols);
for b=1:numel(blocks)
    if ~all(isfield(blocks(b),{'HestimateForDecision','posteriorVarianceByTx'}))
        error('assemble_burst_channel_estimates_v15_1:MissingEstimate', ...
            'Every block must contain its causal channel estimate and posterior variance.');
    end
    idx=blockDataIndices{b};
    if ~isnumeric(idx) || isempty(idx) || ~isvector(idx) || ~isreal(idx) || ...
            any(~isfinite(idx)) || any(idx~=floor(idx)) || any(idx<1) || any(idx>nSymbols) || ...
            any(dataMask(idx))
        error('assemble_burst_channel_estimates_v15_1:InvalidIndices', ...
            'Block data indices must be unique valid symbol indices.');
    end
    localH=blocks(b).HestimateForDecision;
    localPost=blocks(b).posteriorVarianceByTx;
    expectedSize=[nFFT nRx nTx numel(idx)];
    if ~isequal(size(localH),expectedSize) || ...
            ~isequal(size(localPost),[nTx numel(idx)]) || ...
            any(~isfinite(real(localH(:)))) || any(~isfinite(imag(localH(:)))) || ...
            any(~isfinite(localPost(:))) || any(localPost(:)<0)
        error('assemble_burst_channel_estimates_v15_1:EstimateShape', ...
            'Block estimates/posterior variances have invalid dimensions or values.');
    end
    Hest(:,:,:,idx)=localH;
    posteriorVariance(:,idx)=localPost;
    dataMask(idx)=true;
end
end
