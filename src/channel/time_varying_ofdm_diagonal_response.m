function Hdiag = time_varying_ofdm_diagonal_response(tapSequence,nFFT,cpLength,nSymbols)
%TIME_VARYING_OFDM_DIAGONAL_RESPONSE Average per-symbol diagonal OFDM channel.
% For each OFDM symbol, the channel taps are averaged over the useful FFT
% interval (after CP removal), then transformed to frequency domain. This
% produces the diagonal (same-subcarrier) part of the time-varying OFDM
% channel. Off-diagonal inter-carrier coupling is measured separately.
%
% tapSequence: [totalSamples x nTaps x nRx x nTx].
% Hdiag:       [nFFT x nRx x nTx x nSymbols].

if ~isnumeric(tapSequence) || isempty(tapSequence) || ndims(tapSequence) > 4 || ...
        any(~isfinite(real(tapSequence(:)))) || any(~isfinite(imag(tapSequence(:))))
    error('time_varying_ofdm_diagonal_response:InvalidTaps', ...
        'tapSequence must be finite numeric data.');
end
validateInteger(nFFT,'nFFT',2);
validateInteger(cpLength,'cpLength',0);
validateInteger(nSymbols,'nSymbols',1);
blockLength = nFFT+cpLength;
if size(tapSequence,1) ~= blockLength*nSymbols
    error('time_varying_ofdm_diagonal_response:LengthMismatch', ...
        'tapSequence length must equal (nFFT+cpLength)*nSymbols.');
end
nTaps = size(tapSequence,2);
nRx = size(tapSequence,3);
nTx = size(tapSequence,4);
if cpLength < nTaps-1
    error('time_varying_ofdm_diagonal_response:CPTooShort', ...
        'Cyclic prefix must cover the channel memory.');
end

Hdiag = complex(zeros(nFFT,nRx,nTx,nSymbols));
for s = 1:nSymbols
    firstUseful = (s-1)*blockLength+cpLength+1;
    usefulIndices = firstUseful:firstUseful+nFFT-1;
    averageTaps = reshape(mean(tapSequence(usefulIndices,:,:,:),1),nTaps,nRx,nTx);
    Hdiag(:,:,:,s) = fft(averageTaps,nFFT,1);
end
end

function validateInteger(value,name,minValue)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < minValue || value ~= floor(value)
    error('time_varying_ofdm_diagonal_response:InvalidInput', ...
        '%s must be an integer >= %d.',name,minValue);
end
end
