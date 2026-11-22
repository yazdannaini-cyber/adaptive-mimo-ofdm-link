function Hblock = build_time_varying_ofdm_ici_matrix(tapSequence,nFFT,cpLength,symbolIndex)
%BUILD_TIME_VARYING_OFDM_ICI_MATRIX Exact OFDM frequency-coupling matrix.
% Hblock maps all transmit-subcarrier symbols in one OFDM symbol to all
% receive-subcarrier observations:
%   yVec = Hblock*xVec
% where xVec=reshape(Xfd,[nFFT*nTx,1]) and
%       yVec=reshape(Yfd,[nFFT*nRx,1]).
%
% The matrix includes both the diagonal channel and off-diagonal ICI terms
% caused by channel variation within the useful FFT interval. The CP must
% cover the FIR memory so samples before the useful interval wrap to the
% cyclic extension of the same OFDM symbol.
%
% tapSequence dimensions: [totalSamples x nTaps x nRx x nTx].
% Hblock dimensions:      [nFFT*nRx x nFFT*nTx].
%
% This function returns an oracle channel matrix from the true time-varying
% taps. In V14 it is used as a diagnostic/upper-bound receiver reference;
% the mapped pilot pattern does not yet estimate this full matrix.

if ~isnumeric(tapSequence) || isempty(tapSequence) || ndims(tapSequence)>4 || ...
        any(~isfinite(real(tapSequence(:)))) || any(~isfinite(imag(tapSequence(:))))
    error('build_time_varying_ofdm_ici_matrix:InvalidTaps', ...
        'tapSequence must be a finite numeric 4-D array.');
end
validateInteger(nFFT,'nFFT',2);
validateInteger(cpLength,'cpLength',0);
validateInteger(symbolIndex,'symbolIndex',1);
[nTotalSamples,nTaps,nRx,nTx] = size(tapSequence);
blockLength = nFFT+cpLength;
if cpLength < nTaps-1
    error('build_time_varying_ofdm_ici_matrix:CPTooShort', ...
        'cpLength must be at least nTaps-1.');
end
if nTotalSamples < blockLength*symbolIndex
    error('build_time_varying_ofdm_ici_matrix:SymbolOutOfRange', ...
        'tapSequence does not contain the requested OFDM symbol.');
end

firstUseful = (symbolIndex-1)*blockLength+cpLength+1;
usefulIndices = firstUseful:firstUseful+nFFT-1;
Hblock = complex(zeros(nFFT*nRx,nFFT*nTx));

% For output subcarrier k and input subcarrier q, the exact coupling is
%   H(k,q) = (1/N) sum_l exp(-j*2*pi*q*l/N) *
%                    sum_m h_l[m] exp(-j*2*pi*(k-q)*m/N).
% The inner sum is an N-point FFT of the tap trajectory over the useful
% interval. Vectorizing the circulant index avoids a slow four-deep loop.
[kGrid,qGrid] = ndgrid(0:nFFT-1,0:nFFT-1);
fftIndex = mod(kGrid-qGrid,nFFT)+1;
qPhaseIndex = (0:nFFT-1).';
for rx = 1:nRx
    rowIndices = (rx-1)*nFFT+(1:nFFT);
    for tx = 1:nTx
        colIndices = (tx-1)*nFFT+(1:nFFT);
        tapBlock = reshape(tapSequence(usefulIndices,:,rx,tx),nFFT,nTaps);
        coupling = complex(zeros(nFFT,nFFT));
        for ell = 1:nTaps
            tapFFT = fft(tapBlock(:,ell),nFFT);
            delayPhase = exp(-1i*2*pi*qPhaseIndex*(ell-1)/nFFT).';
            coupling = coupling + bsxfun(@times,tapFFT(fftIndex),delayPhase)/nFFT;
        end
        Hblock(rowIndices,colIndices) = coupling;
    end
end
end

function validateInteger(value,name,minValue)
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < minValue || value ~= floor(value)
    error('build_time_varying_ofdm_ici_matrix:InvalidInput', ...
        '%s must be an integer >= %d.',name,minValue);
end
end
