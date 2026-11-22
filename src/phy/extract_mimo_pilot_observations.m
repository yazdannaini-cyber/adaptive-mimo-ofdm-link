function observationsByTx = extract_mimo_pilot_observations(Yfd,pilotSymbols,schedule)
%EXTRACT_MIMO_PILOT_OBSERVATIONS Extract per-link LS observations from mapped pilots.
% Yfd dimensions:       [nFFT x nRx x nSymbols].
% pilotSymbols:         [nFFT x nTx x nPilotOccasions].
% observationsByTx{t}:  [nFFT x nRx x nPilotOccasions].
%
% Since only one transmitter is active in each dedicated pilot symbol,
% element-wise division by the known pilot sequence yields a full-band LS
% observation for all receive antennas and that transmit antenna.

if ~isnumeric(Yfd) || isempty(Yfd) || ndims(Yfd) > 3 || ...
        any(~isfinite(real(Yfd(:)))) || any(~isfinite(imag(Yfd(:))))
    error('extract_mimo_pilot_observations:InvalidYfd', ...
        'Yfd must be a finite numeric [nFFT x nRx x nSymbols] array.');
end
if ~isnumeric(pilotSymbols) || isempty(pilotSymbols) || ndims(pilotSymbols) > 3 || ...
        any(~isfinite(real(pilotSymbols(:)))) || any(~isfinite(imag(pilotSymbols(:))))
    error('extract_mimo_pilot_observations:InvalidPilots', ...
        'pilotSymbols must be a finite numeric [nFFT x nTx x nPilotOccasions] array.');
end
required = {'nSymbols','nTx','pilotIndicesByTx','numPilotSymbols'};
if ~isstruct(schedule) || ~all(isfield(schedule,required)) || ...
        size(Yfd,3) ~= schedule.nSymbols || size(pilotSymbols,2) ~= schedule.nTx
    error('extract_mimo_pilot_observations:InvalidSchedule', ...
        'Schedule dimensions do not match Yfd and pilotSymbols.');
end
if size(Yfd,1) ~= size(pilotSymbols,1)
    error('extract_mimo_pilot_observations:SubcarrierMismatch', ...
        'Yfd and pilotSymbols must have the same number of subcarriers.');
end

nTx = schedule.nTx;
observationsByTx = cell(nTx,1);
for tx = 1:nTx
    pilotIndices = schedule.pilotIndicesByTx{tx};
    if numel(pilotIndices) ~= size(pilotSymbols,3)
        error('extract_mimo_pilot_observations:PilotCountMismatch', ...
            'Each transmitter must have one known pilot sequence per pilot occasion.');
    end
    observations = complex(zeros(size(Yfd,1),size(Yfd,2),numel(pilotIndices)));
    for q = 1:numel(pilotIndices)
        reference = pilotSymbols(:,tx,q);
        if any(abs(reference) <= eps)
            error('extract_mimo_pilot_observations:ZeroPilot', ...
                'Pilot references must be nonzero on every subcarrier.');
        end
        observations(:,:,q) = Yfd(:,:,pilotIndices(q))./reference;
    end
    observationsByTx{tx} = observations;
end
end
