function setup = temporal_channel_tracker_setup(pilotIndices, nSymbols, ...
    pilotNoiseVariance, maxDopplerHz, symbolPeriodSeconds, diffuseVariance)
%TEMPORAL_CHANNEL_TRACKER_SETUP Precompute causal Jakes-LMMSE prediction weights.
% The predictor uses only pilot observations at or before each target symbol.
% The diffuse temporal covariance is R_h(dt)=diffuseVariance*J0(2*pi*fD*dt).
%
% Inputs:
%   pilotIndices          1-based OFDM-symbol indices of pilot observations
%   nSymbols              total number of OFDM symbols in the frame
%   pilotNoiseVariance    complex LS-observation noise variance
%   maxDopplerHz          maximum diffuse Doppler frequency
%   symbolPeriodSeconds   time between adjacent OFDM-symbol snapshots
%   diffuseVariance      per-frequency-coefficient diffuse channel variance
%
% The returned setup is independent of the actual pilot values and can be
% reused across many paired Monte-Carlo frames with the same configuration.

if ~isnumeric(pilotIndices) || isempty(pilotIndices) || ~isvector(pilotIndices) || ~isreal(pilotIndices) || ...
        any(~isfinite(pilotIndices(:))) || any(pilotIndices(:) ~= floor(pilotIndices(:))) || ...
        any(pilotIndices(:) < 1) || any(diff(pilotIndices(:)) <= 0)
    error('temporal_channel_tracker_setup:InvalidPilotIndices', ...
        'pilotIndices must be strictly increasing positive integer indices.');
end
if ~isscalar(nSymbols) || ~isreal(nSymbols) || ~isfinite(nSymbols) || ...
        nSymbols < 1 || nSymbols ~= floor(nSymbols) || pilotIndices(end) > nSymbols
    error('temporal_channel_tracker_setup:InvalidSymbolCount', ...
        'nSymbols must be a positive integer covering all pilot indices.');
end
if ~isnumeric(pilotNoiseVariance) || ~isreal(pilotNoiseVariance) || ...
        isempty(pilotNoiseVariance) || any(~isfinite(pilotNoiseVariance(:))) || ...
        any(pilotNoiseVariance(:)<0) || ...
        (~isscalar(pilotNoiseVariance) && numel(pilotNoiseVariance)~=numel(pilotIndices))
    error('temporal_channel_tracker_setup:InvalidPilotNoiseVariance', ...
        'pilotNoiseVariance must be a nonnegative scalar or one value per pilot.');
end
if isscalar(pilotNoiseVariance)
    noiseByPilot = repmat(double(pilotNoiseVariance),1,numel(pilotIndices));
else
    noiseByPilot = double(pilotNoiseVariance(:).');
end
validateNonnegativeScalar(maxDopplerHz,'maxDopplerHz');
validateNonnegativeScalar(diffuseVariance,'diffuseVariance');
if ~isscalar(symbolPeriodSeconds) || ~isreal(symbolPeriodSeconds) || ...
        ~isfinite(symbolPeriodSeconds) || symbolPeriodSeconds <= 0
    error('temporal_channel_tracker_setup:InvalidSymbolPeriod', ...
        'symbolPeriodSeconds must be finite and positive.');
end

pilotIndices = double(pilotIndices(:).');
setup.nSymbols = double(nSymbols);
setup.pilotIndices = pilotIndices;
setup.pilotNoiseVariance = double(pilotNoiseVariance);
setup.pilotNoiseVarianceByPilot = noiseByPilot;
setup.maxDopplerHz = double(maxDopplerHz);
setup.symbolPeriodSeconds = double(symbolPeriodSeconds);
setup.diffuseVariance = double(diffuseVariance);
setup.observationPositions = cell(nSymbols,1);
setup.weights = cell(nSymbols,1);
setup.posteriorVariance = zeros(nSymbols,1);

for target = 1:nSymbols
    available = find(pilotIndices <= target);
    setup.observationPositions{target} = available;
    if isempty(available) || diffuseVariance == 0
        setup.weights{target} = zeros(numel(available),1);
        if isempty(available)
            setup.posteriorVariance(target) = diffuseVariance;
        else
            setup.posteriorVariance(target) = 0;
        end
        continue;
    end

    observedTimes = pilotIndices(available);
    deltaPP = abs(observedTimes(:)-observedTimes(:).')*symbolPeriodSeconds;
    deltaTP = abs(target-observedTimes(:).')*symbolPeriodSeconds;
    Rpp = diffuseVariance*besselj(0,2*pi*maxDopplerHz*deltaPP);
    rtp = diffuseVariance*besselj(0,2*pi*maxDopplerHz*deltaTP);
    Rpp = (Rpp+Rpp')/2;
    % A tiny diagonal loading handles the rank-one zero-Doppler/noiseless limit.
    loading = max([1,diffuseVariance,noiseByPilot])*1e-12;
    Ryy = Rpp+diag(noiseByPilot(available)+loading);
    weight = Ryy\rtp(:);
    setup.weights{target} = weight;
    setup.posteriorVariance(target) = max(real(diffuseVariance-rtp*weight),0);
end
end

function validateNonnegativeScalar(value,name)
if ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value < 0
    error('temporal_channel_tracker_setup:InvalidScalar', ...
        '%s must be a finite nonnegative scalar.',name);
end
end
