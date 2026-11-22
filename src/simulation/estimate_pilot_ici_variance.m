function [iciNoiseVariance, detail] = estimate_pilot_ici_variance( ...
    observationsByTx, schedule, blockIndex, priorMean, thermalNoiseVariance, ...
    pilotPower, maxDopplerHz, symbolPeriodSeconds, diffuseVariance, ...
    subcarrierSpacingHz, channelEstimateForFloor, previousIciNoiseByBlock)
%ESTIMATE_PILOT_ICI_VARIANCE Causal pilot-innovation estimate of excess ICI.
% For every Tx pilot at the current occasion, form a one-step prediction
% using only earlier pilots. The innovation variance is compared with the
% Jakes prediction-error variance plus known LS thermal-noise variance.
% Positive excess is attributed to unmodelled pilot contamination/ICI.
%
% At the first pilot occasion no temporal innovation is available. In that
% case a small-Doppler analytical ICI floor is used:
%   ICI fraction ~= min(0.1, (pi^2/3)*(fD/DeltaF)^2).
% This is a receiver-side estimate; no true channel or tap sequence is used.

if ~iscell(observationsByTx) || isempty(observationsByTx) || ...
        ~isstruct(schedule) || ~isfield(schedule,'pilotIndicesByTx')
    error('estimate_pilot_ici_variance:InvalidInput', ...
        'Pilot observations and schedule are required.');
end
nTx = numel(observationsByTx);
if blockIndex < 1 || blockIndex ~= floor(blockIndex) || ...
        blockIndex > numel(schedule.pilotOccasionStarts)
    error('estimate_pilot_ici_variance:InvalidBlock','Invalid pilot occasion index.');
end
validatePositive(thermalNoiseVariance,'thermalNoiseVariance');
validatePositive(pilotPower,'pilotPower');
validateNonnegative(maxDopplerHz,'maxDopplerHz');
validatePositive(symbolPeriodSeconds,'symbolPeriodSeconds');
validateNonnegative(diffuseVariance,'diffuseVariance');
validatePositive(subcarrierSpacingHz,'subcarrierSpacingHz');

if nargin < 12 || isempty(previousIciNoiseByBlock)
    previousIciNoiseByBlock = zeros(1,max(blockIndex-1,0));
end
previousIciNoiseByBlock = double(previousIciNoiseByBlock(:).');
if numel(previousIciNoiseByBlock) ~= blockIndex-1 || ...
        any(~isfinite(previousIciNoiseByBlock)) || any(previousIciNoiseByBlock<0)
    error('estimate_pilot_ici_variance:PreviousICIHistory', ...
        'previousIciNoiseByBlock must contain one nonnegative value for each earlier block.');
end
pilotLSNoiseVariance = pilot_ls_noise_profile(thermalNoiseVariance, ...
    pilotPower,previousIciNoiseByBlock);
innovationObserved = nan(nTx,1);
innovationExpected = nan(nTx,1);
excessByTx = zeros(nTx,1);
for tx = 1:nTx
    pilotIndices = schedule.pilotIndicesByTx{tx};
    currentIndex = pilotIndices(blockIndex);
    if blockIndex <= 1
        continue;
    end
    previousIndices = pilotIndices(1:blockIndex-1);
    obsAll = observationsByTx{tx};
    obsPrevious = reshape(obsAll(:,:,1:blockIndex-1), ...
        size(obsAll,1),size(obsAll,2),1,blockIndex-1);

    % observationsByTx{tx} contains one isolated transmit link, whereas
    % priorMean may be either a common single-link prior [N x Nr x 1 x Ns]
    % (legacy/test use) or the full MIMO LOS tensor [N x Nr x Nt x Ns].
    % Select the matching Tx slice for a full MIMO prior. Passing the
    % unsliced full tensor creates an Nt-versus-1 mismatch in the tracker.
    priorTxCount = size(priorMean,3);
    if ~isnumeric(priorMean) || ndims(priorMean)>4 || ...
            size(priorMean,1)~=size(obsAll,1) || ...
            size(priorMean,2)~=size(obsAll,2) || ...
            ~(priorTxCount==1 || priorTxCount==nTx) || ...
            size(priorMean,4)<pilotIndices(end) || ...
            any(~isfinite(real(priorMean(:)))) || ...
            any(~isfinite(imag(priorMean(:))))
        error('estimate_pilot_ici_variance:PriorMeanDimensions', ...
            ['priorMean must be finite [nFFT x nRx x 1 x nSymbols] or ', ...
             '[nFFT x nRx x nTx x nSymbols] data matching the pilot links.']);
    end
    if priorTxCount==1
        priorTxIndex=1;
    else
        priorTxIndex=tx;
    end
    priorMeanTx = reshape(priorMean(:,:,priorTxIndex,:), ...
        size(priorMean,1),size(priorMean,2),1,size(priorMean,4));
    nSymbols = size(priorMeanTx,4);
    setup = temporal_channel_tracker_setup(previousIndices, ...
        nSymbols,pilotLSNoiseVariance,maxDopplerHz, ...
        symbolPeriodSeconds,diffuseVariance);
    [predictedTrack,predictedVariance] = temporal_channel_tracking_apply( ...
        obsPrevious,previousIndices,nSymbols,'jakes-lmmse',priorMeanTx,setup);
    innovation = obsAll(:,:,blockIndex)-predictedTrack(:,:,1,currentIndex);
    observed = mean(abs(innovation(:)).^2);
    expected = thermalNoiseVariance/pilotPower+predictedVariance(currentIndex);
    innovationObserved(tx) = observed;
    innovationExpected(tx) = expected;
    excessByTx(tx) = max(observed-expected,0);
end

channelPower = mean(abs(channelEstimateForFloor(:)).^2);
if ndims(channelEstimateForFloor) >= 3
    summedChannelPower = sum(abs(channelEstimateForFloor).^2,3);
    channelPower = mean(summedChannelPower(:));
end
normalizedDoppler = maxDopplerHz/subcarrierSpacingHz;
iciFractionFloor = min(0.1,(pi^2/3)*normalizedDoppler^2);
floorVariance = iciFractionFloor*max(channelPower,eps);
if blockIndex == 1
    iciNoiseVariance = floorVariance;
else
    iciNoiseVariance = max(sum(excessByTx),floorVariance);
end

if ~isfinite(iciNoiseVariance) || iciNoiseVariance < 0
    error('estimate_pilot_ici_variance:NumericalFailure','Estimated ICI variance is invalid.');
end
detail.innovationObservedVarianceByTx = innovationObserved;
detail.innovationExpectedVarianceByTx = innovationExpected;
detail.excessVarianceByTx = excessByTx;
detail.analyticalICIFloor = floorVariance;
detail.analyticalICIFraction = iciFractionFloor;
detail.estimatedICINoiseVariance = iciNoiseVariance;
end

function validatePositive(x,name)
if ~isscalar(x) || ~isreal(x) || ~isfinite(x) || x <= 0
    error('estimate_pilot_ici_variance:InvalidScalar','%s must be positive and finite.',name);
end
end
function validateNonnegative(x,name)
if ~isscalar(x) || ~isreal(x) || ~isfinite(x) || x < 0
    error('estimate_pilot_ici_variance:InvalidScalar','%s must be nonnegative and finite.',name);
end
end
