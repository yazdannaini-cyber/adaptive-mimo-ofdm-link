function [Hhat, posteriorVariance] = temporal_channel_tracking_apply( ...
    pilotObservations, pilotIndices, nSymbols, method, priorMean, setup)
%TEMPORAL_CHANNEL_TRACKING_APPLY Estimate channel coefficients over OFDM time.
% pilotObservations dimensions: [nFFT x nRx x nTx x nPilotObservations].
% priorMean dimensions:         [nFFT x nRx x nTx x nSymbols], or [].
% Hhat dimensions:              [nFFT x nRx x nTx x nSymbols].
%
% Methods:
%   'hold-last'  causal last-pilot hold, applied to the diffuse residual.
%   'linear'     linear interpolation/extrapolation of the diffuse residual;
%                this is a noncausal offline reference between pilots.
%   'jakes-lmmse'causal model-based predictor using the Jakes J0 covariance.
%
% A known deterministic prior mean (e.g. configured Rician LOS trajectory) is
% removed before tracking and restored at each target symbol for all methods.
% posteriorVariance is an nSymbols vector for Jakes-LMMSE and NaN otherwise.

if ~isnumeric(pilotObservations) || isempty(pilotObservations) || ndims(pilotObservations) > 4 || ...
        any(~isfinite(real(pilotObservations(:)))) || ...
        any(~isfinite(imag(pilotObservations(:))))
    error('temporal_channel_tracking_apply:InvalidObservations', ...
        'pilotObservations must be a nonempty finite numeric array.');
end
if ~isnumeric(pilotIndices) || isempty(pilotIndices) || ~isvector(pilotIndices) || ~isreal(pilotIndices) || ...
        any(~isfinite(pilotIndices(:))) || any(pilotIndices(:) ~= floor(pilotIndices(:))) || ...
        any(pilotIndices(:) < 1) || any(diff(pilotIndices(:)) <= 0)
    error('temporal_channel_tracking_apply:InvalidPilotIndices', ...
        'pilotIndices must be strictly increasing positive integer indices.');
end
if ~isscalar(nSymbols) || ~isreal(nSymbols) || ~isfinite(nSymbols) || ...
        nSymbols < 1 || nSymbols ~= floor(nSymbols) || pilotIndices(end) > nSymbols
    error('temporal_channel_tracking_apply:InvalidSymbolCount', ...
        'nSymbols must cover all pilot indices and be a positive integer.');
end
if ~(ischar(method) || (isstring(method) && isscalar(method)))
    error('temporal_channel_tracking_apply:InvalidMethod', ...
        'method must be hold-last, linear, or jakes-lmmse.');
end
methodText = char(method);
% Canonicalize method names without relying on LOWER.  The interpolation
% branch uses a local pilot-position index; avoiding a variable named
% 'lower' also prevents MATLAB from shadowing the LOWER function.
if strcmpi(methodText,'hold-last')
    method = 'hold-last';
elseif strcmpi(methodText,'linear')
    method = 'linear';
elseif strcmpi(methodText,'jakes-lmmse')
    method = 'jakes-lmmse';
else
    error('temporal_channel_tracking_apply:UnknownMethod', ...
        'Unknown method. Use hold-last, linear, or jakes-lmmse.');
end

[nFFT,nRx,nTx,nPilots] = size(pilotObservations);
pilotIndices = double(pilotIndices(:).');
if nPilots ~= numel(pilotIndices)
    error('temporal_channel_tracking_apply:PilotDimensionMismatch', ...
        'The fourth observation dimension must match numel(pilotIndices).');
end
if isempty(priorMean)
    priorMean = complex(zeros(nFFT,nRx,nTx,nSymbols));
elseif ~isnumeric(priorMean) || ndims(priorMean) > 4 || ...
        any(~isfinite(real(priorMean(:)))) || any(~isfinite(imag(priorMean(:)))) || ...
        size(priorMean,1) ~= nFFT || size(priorMean,2) ~= nRx || ...
        size(priorMean,3) ~= nTx || size(priorMean,4) ~= nSymbols
    error('temporal_channel_tracking_apply:InvalidPriorMean', ...
        'priorMean must be empty or a finite [nFFT x nRx x nTx x nSymbols] array.');
end

pilotMean = priorMean(:,:,:,pilotIndices);
residualPilots = pilotObservations-pilotMean;
Hhat = complex(zeros(nFFT,nRx,nTx,nSymbols));
posteriorVariance = nan(nSymbols,1);

switch method
    case 'hold-last'
        for target = 1:nSymbols
            available = find(pilotIndices <= target,1,'last');
            if isempty(available)
                available = 1; % startup fallback; normal frames place a pilot at symbol 1.
            end
            Hhat(:,:,:,target) = priorMean(:,:,:,target)+residualPilots(:,:,:,available);
        end

    case 'linear'
        for target = 1:nSymbols
            if target <= pilotIndices(1) || numel(pilotIndices) == 1
                residualEstimate = residualPilots(:,:,:,1);
            elseif target >= pilotIndices(end)
                if numel(pilotIndices) >= 2
                    dt = pilotIndices(end)-pilotIndices(end-1);
                    slope = (residualPilots(:,:,:,end)-residualPilots(:,:,:,end-1))/dt;
                    residualEstimate = residualPilots(:,:,:,end) + ...
                        (target-pilotIndices(end))*slope;
                else
                    residualEstimate = residualPilots(:,:,:,end);
                end
            else
                lowerPilotPos = find(pilotIndices <= target,1,'last');
                upperPilotPos = lowerPilotPos+1;
                alpha = (target-pilotIndices(lowerPilotPos))/ ...
                    (pilotIndices(upperPilotPos)-pilotIndices(lowerPilotPos));
                residualEstimate = (1-alpha)*residualPilots(:,:,:,lowerPilotPos) + ...
                    alpha*residualPilots(:,:,:,upperPilotPos);
            end
            Hhat(:,:,:,target) = priorMean(:,:,:,target)+residualEstimate;
        end

    case 'jakes-lmmse'
        if nargin < 6 || ~isstruct(setup) || ...
                ~all(isfield(setup,{'weights','nSymbols','pilotIndices', ...
                'observationPositions','posteriorVariance','diffuseVariance'})) || ...
                setup.nSymbols ~= nSymbols || ~isequal(setup.pilotIndices,pilotIndices)
            error('temporal_channel_tracking_apply:InvalidSetup', ...
                'A matching temporal_channel_tracker_setup result is required.');
        end
        for target = 1:nSymbols
            available = setup.observationPositions{target};
            if isempty(available)
                residualEstimate = complex(zeros(nFFT,nRx,nTx));
            elseif setup.diffuseVariance == 0
                residualEstimate = complex(zeros(nFFT,nRx,nTx));
            else
                residualMatrix = reshape(residualPilots(:,:,:,available), ...
                    nFFT*nRx*nTx,numel(available));
                residualEstimate = reshape(residualMatrix*setup.weights{target}, ...
                    nFFT,nRx,nTx);
            end
            Hhat(:,:,:,target) = priorMean(:,:,:,target)+residualEstimate;
            posteriorVariance(target) = setup.posteriorVariance(target);
        end
end
end
