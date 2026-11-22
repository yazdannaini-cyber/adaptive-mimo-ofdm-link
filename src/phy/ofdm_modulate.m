function txWaveform = ofdm_modulate(Xfd, cpLength)
%OFDM_MODULATE IFFT and cyclic-prefix insertion.
% Xfd dimensions: [nFFT x nTx x nOFDMsymbols].

[nFFT, nTx, nSymbols] = size(Xfd);
if cpLength < 0 || cpLength >= nFFT || cpLength ~= floor(cpLength)
    error('ofdm_modulate:InvalidCP', 'CP length must be an integer in [0, nFFT-1].');
end

txWaveform = complex(zeros(nFFT+cpLength, nTx, nSymbols));
for s = 1:nSymbols
    for t = 1:nTx
        x = ifft(Xfd(:,t,s), nFFT);
        if cpLength > 0
            x = [x(end-cpLength+1:end); x]; %#ok<AGROW>
        end
        txWaveform(:,t,s) = x;
    end
end
end
