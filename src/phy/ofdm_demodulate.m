function Yfd = ofdm_demodulate(rxWaveform, nFFT, cpLength)
%OFDM_DEMODULATE Cyclic-prefix removal and FFT.
% rxWaveform dimensions: [(nFFT+cpLength) x nRx x nOFDMsymbols].

[nSamples, nRx, nSymbols] = size(rxWaveform);
if nSamples ~= nFFT + cpLength
    error('ofdm_demodulate:LengthMismatch', 'Input block length must equal nFFT + cpLength.');
end
Yfd = complex(zeros(nFFT, nRx, nSymbols));
for s = 1:nSymbols
    for r = 1:nRx
        useful = rxWaveform(cpLength+1:cpLength+nFFT, r, s);
        Yfd(:,r,s) = fft(useful, nFFT);
    end
end
end
