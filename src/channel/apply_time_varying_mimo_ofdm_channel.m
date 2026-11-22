function rxWaveform = apply_time_varying_mimo_ofdm_channel(txWaveform, tapSequence, noiseVarianceTime)
%APPLY_TIME_VARYING_MIMO_OFDM_CHANNEL Apply a sample-varying MIMO FIR channel.
% V15.2 implementation vectorizes the sample dimension while preserving
% y_r[n] = sum_{t,l} h_{r,t,l}[n] x_t[n-l+1] exactly.
% txWaveform: [samplesPerOFDMSymbol x nTx x nOFDMSymbols].
% tapSequence: [totalSamples x nTaps x nRx x nTx].
% rxWaveform: [samplesPerOFDMSymbol x nRx x nOFDMSymbols].
if ~isnumeric(txWaveform)||isempty(txWaveform)||ndims(txWaveform)>3|| ...
 any(~isfinite(real(txWaveform(:))))||any(~isfinite(imag(txWaveform(:))))
 error('apply_time_varying_mimo_ofdm_channel:InvalidTransmitSignal', ...
  'txWaveform must be a finite numeric 3-D array.');
end
if ~isnumeric(tapSequence)||isempty(tapSequence)||ndims(tapSequence)>4|| ...
 any(~isfinite(real(tapSequence(:))))||any(~isfinite(imag(tapSequence(:))))
 error('apply_time_varying_mimo_ofdm_channel:InvalidTaps','tapSequence must be a finite numeric 4-D array.');
end
if ~isscalar(noiseVarianceTime)||~isreal(noiseVarianceTime)|| ...
 ~isfinite(noiseVarianceTime)||noiseVarianceTime<0
 error('apply_time_varying_mimo_ofdm_channel:InvalidNoise','noiseVarianceTime must be finite and nonnegative.');
end
[nSamplesPerSymbol,nTx,nSymbols]=size(txWaveform);
[nTotalSamples,nTaps,nRx,nTxChannel]=size(tapSequence);
if nTx~=nTxChannel||nTotalSamples~=nSamplesPerSymbol*nSymbols||nTaps<1||nRx<1
 error('apply_time_varying_mimo_ofdm_channel:DimensionMismatch', ...
  'Waveform and time-varying channel dimensions are inconsistent.');
end
xSerial=complex(zeros(nTotalSamples,nTx));
for tx=1:nTx
 xSerial(:,tx)=reshape(txWaveform(:,tx,:),nTotalSamples,1);
end
ySerial=complex(zeros(nTotalSamples,nRx));
% Vectorized over samples: the channel coefficient is still evaluated at the
% output sample n, exactly as in the original scalar-loop implementation.
for ell=1:nTaps
    nValid=nTotalSamples-ell+1;
    if nValid<=0, continue; end
    for rx=1:nRx
        contribution=complex(zeros(nValid,1));
        for tx=1:nTx
            h=reshape(tapSequence(ell:end,ell,rx,tx),nValid,1);
            contribution=contribution+h.*xSerial(1:nValid,tx);
        end
        ySerial(ell:end,rx)=ySerial(ell:end,rx)+contribution;
    end
end
if noiseVarianceTime>0
 noise=sqrt(noiseVarianceTime/2)*(randn(nTotalSamples,nRx)+1i*randn(nTotalSamples,nRx));
 ySerial=ySerial+noise;
end
rxWaveform=complex(zeros(nSamplesPerSymbol,nRx,nSymbols));
for rx=1:nRx
 rxWaveform(:,rx,:)=reshape(ySerial(:,rx),nSamplesPerSymbol,1,nSymbols);
end
end
