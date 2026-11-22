function [Xhat,postSinr] = mimo_ici_equalize(Yfd,Hblock,noiseVariance,method)
%MIMO_ICI_EQUALIZE Full-band linear MIMO-OFDM equalizer with ICI coupling.
% Yfd:       [nFFT x nRx] received FFT-domain OFDM symbol.
% Hblock:    [nFFT*nRx x nFFT*nTx] full frequency-coupling channel matrix.
% noiseVariance: scalar complex-noise variance per receive subcarrier.
% method:    'mmse' or 'zf'.
% Xhat:      [nFFT x nTx] estimated transmitted QAM grid.
%
% Vector ordering is MATLAB column-major:
%   xVec=reshape(Xfd,[nFFT*nTx,1]);
%   yVec=reshape(Yfd,[nFFT*nRx,1]).
%
% A diagonal-gain correction is applied after linear filtering, matching the
% per-tone receiver convention used by mimo_equalize. This routine is
% computationally more expensive than one-tap equalization because it solves
% a full nFFT*nTx dimensional system.

if ~isnumeric(Yfd) || isempty(Yfd) || ndims(Yfd)>2 || ...
        any(~isfinite(real(Yfd(:)))) || any(~isfinite(imag(Yfd(:))))
    error('mimo_ici_equalize:InvalidReceivedSignal', ...
        'Yfd must be a finite numeric [nFFT x nRx] matrix.');
end
if ~isnumeric(Hblock) || isempty(Hblock) || ndims(Hblock)>2 || ...
        any(~isfinite(real(Hblock(:)))) || any(~isfinite(imag(Hblock(:))))
    error('mimo_ici_equalize:InvalidChannel', ...
        'Hblock must be a finite numeric matrix.');
end
[nFFT,nRx] = size(Yfd);
[nRows,nCols] = size(Hblock);
if nRows ~= nFFT*nRx || mod(nCols,nFFT) ~= 0
    error('mimo_ici_equalize:DimensionMismatch', ...
        'Hblock dimensions must be [nFFT*nRx x nFFT*nTx].');
end
nTx = nCols/nFFT;
if ~isscalar(noiseVariance) || ~isreal(noiseVariance) || ...
        ~isfinite(noiseVariance) || noiseVariance < 0
    error('mimo_ici_equalize:InvalidNoise', ...
        'noiseVariance must be a finite nonnegative scalar.');
end
if ~(ischar(method) || (isstring(method) && isscalar(method)))
    error('mimo_ici_equalize:InvalidMethod','method must be ''mmse'' or ''zf''.');
end
methodText = char(method);
if strcmpi(methodText,'mmse')
    method = 'mmse';
elseif strcmpi(methodText,'zf')
    method = 'zf';
else
    error('mimo_ici_equalize:UnknownMethod','method must be ''mmse'' or ''zf''.');
end

xObservation = reshape(Yfd,nFFT*nRx,1);
if strcmp(method,'zf')
    W = pinv(Hblock);
else
    regularizedGram = Hblock'*Hblock + double(noiseVariance)*eye(nCols);
    W = regularizedGram\Hblock';
end
G = W*Hblock;
raw = W*xObservation;
gain = diag(G);
smallGain = abs(gain)<1e-12;
gain(smallGain) = 1e-12;
xCorrected = raw./gain;
Xhat = reshape(xCorrected,nFFT,nTx);

postSinr = zeros(nFFT,nTx);
for tx = 1:nTx
    indices = (tx-1)*nFFT+(1:nFFT);
    for k = 1:nFFT
        streamIndex = indices(k);
        desired = abs(G(streamIndex,streamIndex))^2;
        totalCoupledPower = sum(abs(G(streamIndex,:)).^2);
        interference = max(totalCoupledPower-desired,0);
        outputNoise = double(noiseVariance)*sum(abs(W(streamIndex,:)).^2);
        postSinr(k,tx) = desired/max(interference+outputNoise,eps);
    end
end
end
