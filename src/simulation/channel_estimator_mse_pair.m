function result = channel_estimator_mse_pair(cfg, channelName, pilotSNRdB, ...
    numFrames, numBootstrap, confidenceLevel, seed)
%CHANNEL_ESTIMATOR_MSE_PAIR Paired estimator-only NMSE/MSE validation.
% Per-link and joint spatial-frequency LMMSE receive exactly the same channel
% and orthogonal-pilot noise realization in every frame. No adaptive detector
% or mode selection is involved, so estimator quality is measured directly.

if nargin < 7
    error('channel_estimator_mse_pair:NotEnoughInputs', ...
        'Provide cfg, channel, pilot SNR, frames, bootstrap count, confidence, and seed.');
end
channelName = lower(char(channelName));
if ~any(strcmp(channelName,{'rayleigh','rician'})) || ...
        ~isscalar(pilotSNRdB) || ~isreal(pilotSNRdB) || ~isfinite(pilotSNRdB) || ...
        ~isscalar(numFrames) || ~isreal(numFrames) || ~isfinite(numFrames) || ...
        numFrames < 2 || numFrames ~= floor(numFrames) || ...
        ~isscalar(numBootstrap) || ~isreal(numBootstrap) || ~isfinite(numBootstrap) || ...
        numBootstrap < 200 || numBootstrap ~= floor(numBootstrap) || ...
        ~isscalar(confidenceLevel) || ~isreal(confidenceLevel) || ~isfinite(confidenceLevel) || ...
        confidenceLevel <= 0 || confidenceLevel >= 1 || ...
        ~isscalar(seed) || ~isreal(seed) || ~isfinite(seed) || seed < 0 || seed ~= floor(seed)
    error('channel_estimator_mse_pair:InvalidInput','Invalid estimator-validation settings.');
end

N = cfg.system.nFFT;
cp = cfg.system.cpLength;
Nt = cfg.system.nTx;
Nr = cfg.system.nRx;
Ppilot = cfg.system.pilotPower;
numCoefficients = N*Nr*Nt;
noiseVarianceFD = 10^(-pilotSNRdB/10);
noiseVarianceTime = noiseVarianceFD/N;

[meanChannel,covarianceScale] = channel_estimation_prior(cfg,channelName);
[~,~,Rrx,Rtx] = spatial_correlation_factors(cfg);
jointModel = channel_estimation_joint_lmmse_setup(cfg.channel.pdp, ...
    noiseVarianceFD,Ppilot,N,covarianceScale,Rtx,Rrx);
[perLinkWeight,perLinkPosteriorVariance] = prepare_per_link_filter( ...
    cfg.channel.pdp,noiseVarianceFD,Ppilot,N,covarianceScale);

pilotGrid = complex(zeros(N,Nt,Nt));
for tx = 1:Nt
    pilotGrid(:,tx,tx) = sqrt(Ppilot);
end
pilotTx = ofdm_modulate(pilotGrid,cp);

frameNmsePerLink = zeros(numFrames,1);
frameNmseJoint = zeros(numFrames,1);
frameMsePerLink = zeros(numFrames,1);
frameMseJoint = zeros(numFrames,1);
for frame = 1:numFrames
    rng(seed+frame,'twister');
    taps = generate_mimo_channel(cfg,channelName);
    Htrue = fft(taps,N,1);
    pilotRx = apply_mimo_ofdm_channel(pilotTx,taps,noiseVarianceTime);
    pilotY = ofdm_demodulate(pilotRx,N,cp);
    Hls = channel_estimation_ls(pilotY,Ppilot);

    residualLS = reshape(Hls-meanChannel,N,Nr*Nt);
    filteredResidual = perLinkWeight*residualLS;
    HperLink = meanChannel + reshape(filteredResidual,N,Nr,Nt);
    Hjoint = channel_estimation_joint_lmmse_apply(Hls,meanChannel,jointModel);
    errorPerLink = HperLink-Htrue;
    errorJoint = Hjoint-Htrue;
    channelEnergy = max(sum(abs(Htrue(:)).^2),eps);
    frameNmsePerLink(frame) = sum(abs(errorPerLink(:)).^2)/channelEnergy;
    frameNmseJoint(frame) = sum(abs(errorJoint(:)).^2)/channelEnergy;
    frameMsePerLink(frame) = sum(abs(errorPerLink(:)).^2)/numCoefficients;
    frameMseJoint(frame) = sum(abs(errorJoint(:)).^2)/numCoefficients;
end

seedPerLink = seed+700000;
seedJoint = seed+800000;
seedDifference = seed+900000;
ciPerLink = frame_bootstrap_mean_ci(frameNmsePerLink,numBootstrap,confidenceLevel,seedPerLink);
ciJoint = frame_bootstrap_mean_ci(frameNmseJoint,numBootstrap,confidenceLevel,seedJoint);
pairedMseDifference = frameMsePerLink-frameMseJoint;
ciMseDifference = frame_bootstrap_mean_ci(pairedMseDifference,numBootstrap,confidenceLevel,seedDifference);

nmsePerLink = mean(frameNmsePerLink);
nmseJoint = mean(frameNmseJoint);
msePerLink = mean(frameMsePerLink);
mseJoint = mean(frameMseJoint);
predictedPerLinkMSE = mean(perLinkPosteriorVariance);
predictedJointMSE = mean(jointModel.posteriorMeanErrorVarianceByTone);
if nmsePerLink > 0
    nmseImprovementPercent = 100*(nmsePerLink-nmseJoint)/nmsePerLink;
else
    nmseImprovementPercent = 0;
end
if predictedPerLinkMSE > eps
    perLinkCalibrationRatio = msePerLink/predictedPerLinkMSE;
else
    perLinkCalibrationRatio = NaN;
end
if predictedJointMSE > eps
    jointCalibrationRatio = mseJoint/predictedJointMSE;
else
    jointCalibrationRatio = NaN;
end

rhoTx = cfg.channel.spatialCorrelation.txRho;
rhoRx = cfg.channel.spatialCorrelation.rxRho;
result = table(string(channelName),pilotSNRdB,rhoTx,rhoRx,numFrames, ...
    nmsePerLink,10*log10(max(nmsePerLink,realmin)),ciPerLink(1),ciPerLink(2), ...
    nmseJoint,10*log10(max(nmseJoint,realmin)),ciJoint(1),ciJoint(2), ...
    nmseImprovementPercent,msePerLink,mseJoint,predictedPerLinkMSE, ...
    predictedJointMSE,perLinkCalibrationRatio,jointCalibrationRatio, ...
    mean(pairedMseDifference),ciMseDifference(1),ciMseDifference(2), ...
    'VariableNames',{'Channel','PilotSNR_dB','TxCorrelationCoefficient', ...
    'RxCorrelationCoefficient','NumFrames','PerLink_NMSE','PerLink_NMSE_dB', ...
    'PerLink_NMSE_CI95_Lower','PerLink_NMSE_CI95_Upper','Joint_NMSE', ...
    'Joint_NMSE_dB','Joint_NMSE_CI95_Lower','Joint_NMSE_CI95_Upper', ...
    'Joint_NMSE_Improvement_Percent','PerLink_EmpiricalMSE_PerCoefficient', ...
    'Joint_EmpiricalMSE_PerCoefficient','PerLink_PredictedPosteriorMSE_PerCoefficient', ...
    'Joint_PredictedPosteriorMSE_PerCoefficient','PerLink_MSE_CalibrationRatio', ...
    'Joint_MSE_CalibrationRatio','Paired_MSE_Improvement_PerCoefficient', ...
    'Paired_MSE_Improvement_CI95_Lower','Paired_MSE_Improvement_CI95_Upper'});
end

function [W,posteriorVariance] = prepare_per_link_filter(pdp,noiseVariance,pilotPower,nFFT,covarianceScale)
pdp = pdp(:)/sum(pdp);
k = (0:nFFT-1).';
R = complex(zeros(nFFT,nFFT));
for tap = 0:numel(pdp)-1
    steering = exp(-1i*2*pi*k*tap/nFFT);
    R = R+pdp(tap+1)*(steering*steering');
end
R = covarianceScale*(R+R')/2;
noiseLS = noiseVariance/pilotPower;
if covarianceScale == 0
    W = zeros(nFFT,nFFT);
    posteriorVariance = zeros(nFFT,1);
elseif noiseLS <= eps
    W = R*pinv(R);
    posteriorVariance = zeros(nFFT,1);
else
    W = R/(R+noiseLS*eye(nFFT));
    Rpost = R-W*R;
    Rpost = (Rpost+Rpost')/2;
    posteriorVariance = max(real(diag(Rpost)),0);
end
end
