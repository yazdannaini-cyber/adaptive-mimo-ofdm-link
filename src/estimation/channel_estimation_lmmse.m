function Hlmmse = channel_estimation_lmmse(Hls, pdp, noiseVariance, pilotPower, meanResponse, covarianceScale)
%CHANNEL_ESTIMATION_LMMSE Model-aware frequency-domain LMMSE channel estimator.
% Supports a nonzero prior mean and a scaled diffuse-channel covariance.
% Optional arguments preserve the original zero-mean Rayleigh-style interface:
%   meanResponse    : [nFFT x nRx x nTx] prior channel mean (default zero)
%   covarianceScale : diffuse covariance multiplier (default one)
%
% For a Rician channel, use the deterministic LOS response as meanResponse
% and 1/(K+1) as covarianceScale. For the deterministic AWGN identity channel,
% use meanResponse=I and covarianceScale=0.

[nFFT,nRx,nTx] = size(Hls);
pdp = pdp(:);
if nargin < 5 || isempty(meanResponse)
    meanResponse = complex(zeros(size(Hls)));
end
if nargin < 6 || isempty(covarianceScale)
    covarianceScale = 1;
end
if ~isscalar(noiseVariance) || isnan(noiseVariance) || isinf(noiseVariance) || noiseVariance < 0 || ...
        ~isscalar(pilotPower) || isnan(pilotPower) || isinf(pilotPower) || pilotPower <= 0 || ...
        any(isnan(pdp)) || any(isinf(pdp)) || any(pdp < 0) || sum(pdp) <= 0 || ...
        ~isscalar(covarianceScale) || (isnan(covarianceScale) || isinf(covarianceScale)) || covarianceScale < 0
    error('channel_estimation_lmmse:InvalidInput', ...
        'Invalid noise, pilot power, PDP, or covariance scale.');
end
if any(isnan(meanResponse(:))) || any(isinf(meanResponse(:)))
    error('channel_estimation_lmmse:MeanFinite','meanResponse must contain finite values.');
end
if ~isequal(size(meanResponse),size(Hls))
    error('channel_estimation_lmmse:MeanSize', ...
        'meanResponse must have the same dimensions as Hls.');
end
pdp = pdp/sum(pdp);

% A deterministic prior (for example, the configured AWGN identity channel)
% has no diffuse uncertainty, so its MMSE estimate is exactly the prior mean.
if covarianceScale == 0
    Hlmmse = meanResponse;
    return;
end

subcarrierIndex = (0:nFFT-1).';
R = complex(zeros(nFFT,nFFT));
for tap = 0:numel(pdp)-1
    steering = exp(-1i*2*pi*subcarrierIndex*tap/nFFT);
    R = R + pdp(tap+1)*(steering*steering');
end
R = covarianceScale*R;
noiseLS = (noiseVariance/pilotPower)*eye(nFFT);
if noiseVariance <= eps
    % The covariance matrix is low rank when the channel has fewer taps than
    % subcarriers. Projection preserves every valid noiseless channel response.
    W = R*pinv(R);
else
    W = R/(R+noiseLS);
end

Hlmmse = complex(zeros(size(Hls)));
for r = 1:nRx
    for t = 1:nTx
        residual = Hls(:,r,t)-meanResponse(:,r,t);
        Hlmmse(:,r,t) = meanResponse(:,r,t)+W*residual;
    end
end
end
