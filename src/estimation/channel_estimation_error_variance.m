function errorVariance = channel_estimation_error_variance(pdp, noiseVariance, pilotPower, nFFT, covarianceScale)
%CHANNEL_ESTIMATION_ERROR_VARIANCE Posterior per-link channel error variance.
% Returns diag(R_post), where
%   R_post = R_h - R_h*(R_h + sigma_LS^2 I)^(-1)*R_h.
% The result is the same for each independent MIMO channel link under the
% configured PDP and is used to account for channel-estimation uncertainty
% in adaptive detection and link-quality prediction.

if nargin < 5 || isempty(covarianceScale)
    covarianceScale = 1;
end
pdp = pdp(:);
if nFFT < 1 || nFFT ~= floor(nFFT) || ~isscalar(pilotPower) || ...
        isnan(pilotPower) || isinf(pilotPower) || pilotPower <= 0 || ...
        ~isscalar(noiseVariance) || isnan(noiseVariance) || isinf(noiseVariance) || noiseVariance < 0 || ...
        isempty(pdp) || any(isnan(pdp)) || any(isinf(pdp)) || any(pdp < 0) || sum(pdp) <= 0 || ...
        ~isscalar(covarianceScale) || (isnan(covarianceScale) || isinf(covarianceScale)) || covarianceScale < 0
    error('channel_estimation_error_variance:InvalidInput','Invalid estimator parameters.');
end
if noiseVariance <= eps || covarianceScale == 0
    errorVariance = zeros(nFFT,1);
    return;
end

pdp = pdp/sum(pdp);
subcarrierIndex = (0:nFFT-1).';
R = complex(zeros(nFFT,nFFT));
for tap = 0:numel(pdp)-1
    steering = exp(-1i*2*pi*subcarrierIndex*tap/nFFT);
    R = R + pdp(tap+1)*(steering*steering');
end
R = covarianceScale*R;
noiseLS = (noiseVariance/pilotPower)*eye(nFFT);
W = R/(R+noiseLS);
Rpost = R-W*R;
Rpost = (Rpost+Rpost')/2; % suppress numerical anti-Hermitian residue
errorVariance = max(real(diag(Rpost)),0);
end
