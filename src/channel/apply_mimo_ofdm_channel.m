function rxWaveform = apply_mimo_ofdm_channel(txWaveform, channelTaps, noiseVarianceTime)
%APPLY_MIMO_OFDM_CHANNEL Apply per-link FIR channels and complex AWGN.
% The channel is static over each OFDM block. CP insertion/removal is handled
% by the OFDM functions, and the caller must ensure CP >= channel memory.

[nSamples, nTx, nSymbols] = size(txWaveform);
[nTaps, nRx, nTxChannel] = size(channelTaps);
if nTx ~= nTxChannel
    error('apply_mimo_ofdm_channel:TxMismatch', 'Tx waveform and channel dimensions differ.');
end
if noiseVarianceTime < 0
    error('apply_mimo_ofdm_channel:NegativeNoise', 'Noise variance cannot be negative.');
end

rxWaveform = complex(zeros(nSamples, nRx, nSymbols));
for s = 1:nSymbols
    for r = 1:nRx
        y = complex(zeros(nSamples,1));
        for t = 1:nTx
            h = channelTaps(:,r,t);
            filtered = conv(txWaveform(:,t,s), h);
            y = y + filtered(1:nSamples);
        end
        if noiseVarianceTime > 0
            noise = sqrt(noiseVarianceTime/2) * ...
                (randn(nSamples,1) + 1i*randn(nSamples,1));
            y = y + noise;
        end
        rxWaveform(:,r,s) = y;
    end
end
end
