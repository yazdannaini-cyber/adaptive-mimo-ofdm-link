function schedule = mimo_ofdm_pilot_schedule(nSymbols,nTx,pilotSpacing)
%MIMO_OFDM_PILOT_SCHEDULE Build orthogonal time-division MIMO pilot slots.
% Each pilot occasion reserves nTx consecutive OFDM symbols. During the
% symbol assigned to transmitter t, only that transmitter sends a known
% full-band pilot. The per-transmitter pilot timelines are therefore
% staggered by t-1 symbols and are handled independently by the tracker.

validateIntegerScalar(nSymbols,'nSymbols',1);
validateIntegerScalar(nTx,'nTx',1);
validateIntegerScalar(pilotSpacing,'pilotSpacing',1);
if pilotSpacing < nTx
    error('mimo_ofdm_pilot_schedule:SpacingTooShort', ...
        'pilotSpacing must be at least nTx to avoid overlapping pilot slots.');
end
if nSymbols < nTx
    error('mimo_ofdm_pilot_schedule:FrameTooShort', ...
        'The frame must contain at least nTx symbols for one orthogonal pilot occasion.');
end

starts = 1:pilotSpacing:(nSymbols-nTx+1);
if isempty(starts)
    error('mimo_ofdm_pilot_schedule:NoPilotOccasion', ...
        'No complete orthogonal pilot occasion fits in the frame.');
end
pilotIndicesByTx = cell(nTx,1);
pilotMask = false(1,nSymbols);
for tx = 1:nTx
    idx = starts+(tx-1);
    pilotIndicesByTx{tx} = idx;
    pilotMask(idx) = true;
end
dataIndices = find(~pilotMask);
if isempty(dataIndices)
    error('mimo_ofdm_pilot_schedule:NoDataSymbols', ...
        'The pilot schedule leaves no data OFDM symbols; increase pilotSpacing or frame length.');
end

schedule.nSymbols = nSymbols;
schedule.nTx = nTx;
schedule.pilotSpacing = pilotSpacing;
schedule.pilotOccasionStarts = starts;
schedule.pilotIndicesByTx = pilotIndicesByTx;
schedule.pilotMask = pilotMask;
schedule.allPilotIndices = find(pilotMask);
schedule.dataIndices = dataIndices;
schedule.numPilotSymbols = nnz(pilotMask);
schedule.numDataSymbols = numel(dataIndices);
schedule.pilotOverhead = schedule.numPilotSymbols/nSymbols;
end

function validateIntegerScalar(value,name,minValue)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < minValue || value ~= floor(value)
    error('mimo_ofdm_pilot_schedule:InvalidInput', ...
        '%s must be an integer greater than or equal to %d.',name,minValue);
end
end
