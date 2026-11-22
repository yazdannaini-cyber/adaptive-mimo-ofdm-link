function covarianceByTone = channel_error_effective_noise_covariance(model, thermalNoiseVariance, streamPower)
%CHANNEL_ERROR_EFFECTIVE_NOISE_COVARIANCE Build colored receive-noise covariance.
% The channel-estimation error is E(:,:,k), and the received error term is
% E(:,:,k)*x. For independent, unit-power transmit streams,
%   R_eff(k) = sigma_n^2 I + streamPower * sum_t Cov(e_t(k)),
% where Cov(e_t) is the receive-vector error covariance for Tx column t.
% Cross-Tx error covariance does not contribute to this second moment when
% the transmitted streams are independent and zero-mean.
%
% Inputs:
%   model.posteriorLinkCovarianceByTone: [nRx*nTx x nRx*nTx x nFFT]
%   thermalNoiseVariance: nonnegative scalar or nFFT-element vector
%   streamPower: nonnegative scalar (default 1)
%
% Output:
%   covarianceByTone: [nRx x nRx x nFFT], Hermitian positive definite when
%   thermalNoiseVariance is positive.

if nargin < 3 || isempty(streamPower)
    streamPower = 1;
end
required = {'nFFT','nRx','nTx','posteriorLinkCovarianceByTone'};
if ~isstruct(model) || ~all(isfield(model,required))
    error('channel_error_effective_noise_covariance:InvalidModel', ...
        'A joint-LMMSE model with full posterior link covariance is required.');
end
if ~isscalar(streamPower) || ~isreal(streamPower) || ~isfinite(streamPower) || streamPower < 0
    error('channel_error_effective_noise_covariance:InvalidStreamPower', ...
        'streamPower must be a finite nonnegative scalar.');
end
if isscalar(thermalNoiseVariance)
    noiseByTone = repmat(double(thermalNoiseVariance),model.nFFT,1);
elseif isnumeric(thermalNoiseVariance) && numel(thermalNoiseVariance) == model.nFFT
    noiseByTone = double(thermalNoiseVariance(:));
else
    error('channel_error_effective_noise_covariance:NoiseSize', ...
        'Thermal noise must be scalar or contain one value per tone.');
end
if ~isreal(noiseByTone) || any(~isfinite(noiseByTone)) || any(noiseByTone < 0)
    error('channel_error_effective_noise_covariance:InvalidNoise', ...
        'Thermal noise variances must be finite and nonnegative.');
end

nLinks = model.nRx*model.nTx;
expectedSize = [nLinks nLinks model.nFFT];
Cpost = model.posteriorLinkCovarianceByTone;
if ~isequal(size(Cpost),expectedSize) || ...
        any(~isfinite(real(Cpost(:)))) || any(~isfinite(imag(Cpost(:))))
    error('channel_error_effective_noise_covariance:CovarianceSize', ...
        'Posterior link covariance has inconsistent dimensions or non-finite values.');
end

covarianceByTone = complex(zeros(model.nRx,model.nRx,model.nFFT));
I = eye(model.nRx);
for k = 1:model.nFFT
    channelErrorCovariance = complex(zeros(model.nRx));
    for tx = 1:model.nTx
        idx = (tx-1)*model.nRx+(1:model.nRx);
        channelErrorCovariance = channelErrorCovariance + Cpost(idx,idx,k);
    end
    channelErrorCovariance = (channelErrorCovariance+channelErrorCovariance')/2;
    R = noiseByTone(k)*I + streamPower*channelErrorCovariance;
    covarianceByTone(:,:,k) = (R+R')/2;
end
end
