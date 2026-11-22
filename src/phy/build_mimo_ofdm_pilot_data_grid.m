function [Xfd,dataBits,pilotSymbols,schedule] = build_mimo_ofdm_pilot_data_grid( ...
    nFFT,nTx,nSymbols,pilotSpacing,modulationOrder,pilotPower,dataSeed,pilotSeed)
%BUILD_MIMO_OFDM_PILOT_DATA_GRID Map orthogonal pilots and QAM data to OFDM REs.
% Pilot scheme: one transmitter active per dedicated OFDM pilot symbol, with
% known constant-modulus QPSK on every subcarrier. Pilot slots are time-
% division orthogonal across transmit antennas. All remaining REs carry data.
%
% Outputs:
%   Xfd          [nFFT x nTx x nSymbols] frequency-domain grid
%   dataBits     [nFFT*log2(M) x nTx x nDataSymbols] transmitted data bits
%   pilotSymbols [nFFT x nTx x nPilotOccasions] known pilot sequences
%   schedule     indices and pilot-overhead metadata

validatePositiveInteger(nFFT,'nFFT');
validatePositiveInteger(nTx,'nTx');
validatePositiveInteger(nSymbols,'nSymbols');
validatePositiveInteger(modulationOrder,'modulationOrder');
if ~ismember(modulationOrder,[4 16 64])
    error('build_mimo_ofdm_pilot_data_grid:InvalidModulation', ...
        'modulationOrder must be 4, 16, or 64.');
end
if ~isscalar(pilotPower) || ~isreal(pilotPower) || ~isfinite(pilotPower) || pilotPower <= 0
    error('build_mimo_ofdm_pilot_data_grid:InvalidPilotPower', ...
        'pilotPower must be finite and positive.');
end
validateSeed(dataSeed,'dataSeed');
validateSeed(pilotSeed,'pilotSeed');
schedule = mimo_ofdm_pilot_schedule(nSymbols,nTx,pilotSpacing);

oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>

nOccasions = numel(schedule.pilotOccasionStarts);
pilotSymbols = complex(zeros(nFFT,nTx,nOccasions));
Xfd = complex(zeros(nFFT,nTx,nSymbols));

% Deterministic known QPSK pilots, generated independently of payload bits.
rng(pilotSeed,'twister');
for tx = 1:nTx
    for q = 1:nOccasions
        pilotBits = randi([0 1],2*nFFT,1);
        sequence = sqrt(pilotPower)*qam_mod(pilotBits,4);
        pilotSymbols(:,tx,q) = sequence;
        pilotIndex = schedule.pilotIndicesByTx{tx}(q);
        Xfd(:,tx,pilotIndex) = sequence;
    end
end

% Payload data on all subcarriers of every non-pilot OFDM symbol.
rng(dataSeed,'twister');
bitsPerQAM = log2(modulationOrder);
dataBits = zeros(nFFT*bitsPerQAM,nTx,schedule.numDataSymbols);
for q = 1:schedule.numDataSymbols
    symbolIndex = schedule.dataIndices(q);
    for tx = 1:nTx
        bits = randi([0 1],nFFT*bitsPerQAM,1);
        dataBits(:,tx,q) = bits;
        Xfd(:,tx,symbolIndex) = qam_mod(bits,modulationOrder);
    end
end
end

function validatePositiveInteger(value,name)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < 1 || value ~= floor(value)
    error('build_mimo_ofdm_pilot_data_grid:InvalidInput', ...
        '%s must be a positive integer.',name);
end
end

function validateSeed(value,name)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < 0 || value ~= floor(value)
    error('build_mimo_ofdm_pilot_data_grid:InvalidSeed', ...
        '%s must be a nonnegative integer.',name);
end
end
