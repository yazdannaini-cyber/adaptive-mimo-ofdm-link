function [Hest, errorVarianceByTone, model] = channel_estimation_joint_lmmse( ...
    Hls, pdp, noiseVariance, pilotPower, meanResponse, covarianceScale, Rtx, Rrx)
%CHANNEL_ESTIMATION_JOINT_LMMSE One-shot joint spatial-frequency LMMSE estimate.
% For repeated Monte Carlo estimation, prepare the model once with
% channel_estimation_joint_lmmse_setup and apply it per frame.

if nargin < 5, meanResponse = []; end
if nargin < 6, covarianceScale = 1; end
if nargin < 7, Rtx = []; end
if nargin < 8, Rrx = []; end
[nFFT,nRx,nTx] = size(Hls);
if isempty(Rtx), Rtx = eye(nTx); end
if isempty(Rrx), Rrx = eye(nRx); end
if isempty(meanResponse), meanResponse = complex(zeros(size(Hls))); end
model = channel_estimation_joint_lmmse_setup(pdp,noiseVariance,pilotPower, ...
    nFFT,covarianceScale,Rtx,Rrx);
Hest = channel_estimation_joint_lmmse_apply(Hls,meanResponse,model);
errorVarianceByTone = model.posteriorErrorVarianceByTone;
end
