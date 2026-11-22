function Hest = channel_estimation_joint_lmmse_apply(Hls, meanResponse, model)
%CHANNEL_ESTIMATION_JOINT_LMMSE_APPLY Apply a prepared joint LMMSE transform.
% Hls and meanResponse have dimensions [nFFT x nRx x nTx].

if ~isnumeric(Hls) || ndims(Hls) > 3 || isempty(Hls) || ...
        any(~isfinite(real(Hls(:)))) || any(~isfinite(imag(Hls(:))))
    error('channel_estimation_joint_lmmse_apply:InvalidHls', ...
        'Hls must be a finite numeric 3-D array.');
end
if ~isstruct(model) || ~all(isfield(model,{'nFFT','nRx','nTx','Ufreq','Ulink','gain','deterministicPrior'}))
    error('channel_estimation_joint_lmmse_apply:InvalidModel','A prepared joint-LMMSE model is required.');
end
if nargin < 2 || isempty(meanResponse)
    meanResponse = complex(zeros(size(Hls)));
end
if size(Hls,1) ~= model.nFFT || size(Hls,2) ~= model.nRx || ...
        size(Hls,3) ~= model.nTx || size(meanResponse,1) ~= model.nFFT || ...
        size(meanResponse,2) ~= model.nRx || size(meanResponse,3) ~= model.nTx || ...
        any(~isfinite(real(meanResponse(:)))) || any(~isfinite(imag(meanResponse(:))))
    error('channel_estimation_joint_lmmse_apply:DimensionMismatch', ...
        'Hls and meanResponse dimensions must match the prepared model.');
end
if model.deterministicPrior
    Hest = meanResponse;
    return;
end

nLinks = model.nRx*model.nTx;
residual = reshape(Hls-meanResponse,model.nFFT,nLinks);
% vec(Ufreq'*E*conj(Ulink)) = (Ulink'*kron Ufreq')*vec(E).
transformed = model.Ufreq' * residual * conj(model.Ulink);
filtered = model.Ufreq * (model.gain.*transformed) * model.Ulink.';
Hest = meanResponse + reshape(filtered,size(Hls));
end
