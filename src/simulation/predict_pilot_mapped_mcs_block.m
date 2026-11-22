function prediction = predict_pilot_mapped_mcs_block(HbySymbol, thermalNoiseVariance, ...
    modulations, targetBER, spectralEfficiency, posteriorVarianceByTx, iciNoiseBySymbol)
%PREDICT_PILOT_MAPPED_MCS_BLOCK Uncertainty/ICI-aware MCS prediction over a data burst.
% HbySymbol: [nFFT x nRx x 2 x nDataSymbols], estimated channels available
%            at the decision instant (future pilots are excluded by caller).
% posteriorVarianceByTx: [2 x nDataSymbols], temporal-LMMSE posterior error
%                        variance for each transmit-link column and symbol.
% iciNoiseBySymbol: scalar or [1 x nDataSymbols] equivalent receive-noise
%                   variance from causal pilot-innovation analysis.
%
% The predictor deliberately treats channel-estimation error and residual ICI
% as additional uncorrelated receive disturbance. It is conservative and is
% still a link-adaptation estimate, not a guarantee of measured BER.

if ~isnumeric(HbySymbol) || isempty(HbySymbol) || ndims(HbySymbol) > 4 || ...
        size(HbySymbol,3) ~= 2 || ...
        any(~isfinite(real(HbySymbol(:)))) || any(~isfinite(imag(HbySymbol(:))))
    error('predict_pilot_mapped_mcs_block:InvalidChannel', ...
        'HbySymbol must be finite [nFFT x nRx x 2 x nDataSymbols].');
end
if ~isscalar(thermalNoiseVariance) || ~isreal(thermalNoiseVariance) || ...
        ~isfinite(thermalNoiseVariance) || thermalNoiseVariance <= 0
    error('predict_pilot_mapped_mcs_block:InvalidNoise', ...
        'thermalNoiseVariance must be finite and positive.');
end
modulations = double(modulations(:).');
if isempty(modulations) || any(diff(modulations) <= 0) || ...
        any(modulations < 4) || ...
        any(abs(log2(modulations)-round(log2(modulations))) > eps) || ...
        any(abs(log2(sqrt(modulations))-round(log2(sqrt(modulations)))) > eps)
    error('predict_pilot_mapped_mcs_block:InvalidModulations', ...
        'modulations must be ascending square power-of-four QAM orders.');
end
if ~isscalar(targetBER) || ~isreal(targetBER) || ~isfinite(targetBER) || ...
        targetBER <= 0 || targetBER >= 0.5
    error('predict_pilot_mapped_mcs_block:InvalidTarget', ...
        'targetBER must lie in (0,0.5).');
end
[nFFT,nRx,nTx,nData] = size(HbySymbol);
if nTx ~= 2 || nData < 2 || mod(nData,2) ~= 0
    error('predict_pilot_mapped_mcs_block:BlockShape', ...
        'Each block must contain an even number of data symbols and two Tx antennas.');
end
if ~isequal(size(posteriorVarianceByTx),[2 nData]) || ...
        any(~isfinite(posteriorVarianceByTx(:))) || any(posteriorVarianceByTx(:) < 0)
    error('predict_pilot_mapped_mcs_block:PosteriorVariance', ...
        'posteriorVarianceByTx must be finite nonnegative [2 x nDataSymbols].');
end
if isscalar(iciNoiseBySymbol)
    iciNoiseBySymbol = repmat(double(iciNoiseBySymbol),1,nData);
else
    iciNoiseBySymbol = double(iciNoiseBySymbol(:).');
end
if numel(iciNoiseBySymbol) ~= nData || ~isreal(iciNoiseBySymbol) || ...
        any(~isfinite(iciNoiseBySymbol)) || any(iciNoiseBySymbol < 0)
    error('predict_pilot_mapped_mcs_block:IciVariance', ...
        'iciNoiseBySymbol must be a nonnegative scalar or one value per data symbol.');
end
if ~isnumeric(spectralEfficiency) || ~isequal(size(spectralEfficiency),[2 numel(modulations)]) || ...
        any(~isfinite(spectralEfficiency(:))) || any(spectralEfficiency(:) <= 0)
    error('predict_pilot_mapped_mcs_block:InvalidEfficiency', ...
        'spectralEfficiency must be a positive 2-by-numel(modulations) matrix.');
end

% Spatial-multiplexing SINR is evaluated independently for every OFDM slot.
smSinr = zeros(nFFT,2,nData);
effectiveNoise = thermalNoiseVariance + iciNoiseBySymbol;
for s = 1:nData
    Hs = reshape(HbySymbol(:,:,:,s),nFFT,nRx,2);
    effectiveNoise(s) = effectiveNoise(s) + sum(posteriorVarianceByTx(:,s));
    [~,localSinr] = mimo_equalize(complex(zeros(nFFT,nRx,1)),Hs, ...
        effectiveNoise(s),'mmse');
    smSinr(:,:,s) = localSinr(:,:,1);
end

% Alamouti uses the actual two channel predictions in each adjacent slot.
alSinr = zeros(nFFT,2,nData/2);
alamoutiNoise = zeros(1,nData/2);
for pair = 1:nData/2
    s1 = 2*pair-1;
    s2 = 2*pair;
    H1 = reshape(HbySymbol(:,:,:,s1),nFFT,nRx,2);
    H2 = reshape(HbySymbol(:,:,:,s2),nFFT,nRx,2);
    sigma1 = thermalNoiseVariance + iciNoiseBySymbol(s1) + sum(posteriorVarianceByTx(:,s1));
    sigma2 = thermalNoiseVariance + iciNoiseBySymbol(s2) + sum(posteriorVarianceByTx(:,s2));
    alSinr(:,:,pair) = alamouti_sinr_two_slot(H1,H2,sigma1,sigma2);
    alamoutiNoise(pair) = mean([sigma1 sigma2]);
end

predictedBER = zeros(2,numel(modulations));
for m = 1:numel(modulations)
    smValues = qam_ber_approx(smSinr(:),modulations(m));
    alValues = qam_ber_approx(alSinr(:),modulations(m));
    predictedBER(1,m) = mean(smValues(:));
    predictedBER(2,m) = mean(alValues(:));
end
choice = select_mimo_link_mode(predictedBER,modulations,targetBER,spectralEfficiency);

prediction.predictedBER = predictedBER;
prediction.spectralEfficiency = spectralEfficiency;
prediction.modulations = modulations;
prediction.postSinrSpatialMultiplexing = smSinr;
prediction.postSinrAlamouti = alSinr;
prediction.effectiveNoiseVarianceBySymbol = effectiveNoise;
prediction.alamoutiEffectiveNoiseVarianceByPair = alamoutiNoise;
prediction.meanPosteriorVariance = mean(posteriorVarianceByTx(:));
prediction.meanIciNoiseVariance = mean(iciNoiseBySymbol);
prediction.choice = choice;
end

function postSinr = alamouti_sinr_two_slot(Hfirst,Hsecond,noiseFirst,noiseSecond)
[nFFT,nRx,~] = size(Hfirst);
postSinr = zeros(nFFT,2);
for k = 1:nFFT
    H1 = reshape(Hfirst(k,:,:),nRx,2);
    H2 = reshape(Hsecond(k,:,:),nRx,2);
    A = [H1(:,1), H1(:,2); conj(H2(:,2)), -conj(H2(:,1))];
    R = diag([repmat(max(noiseFirst,1e-12),1,nRx), ...
              repmat(max(noiseSecond,1e-12),1,nRx)]);
    RinvA = R\A;
    Rinv = R\eye(2*nRx);
    W = (A'*RinvA+eye(2))\(A'*Rinv);
    G = W*A;
    for stream = 1:2
        desired = abs(G(stream,stream))^2;
        interference = max(sum(abs(G(stream,:)).^2)-desired,0);
        outputNoise = max(real(W(stream,:)*R*W(stream,:)'),0);
        postSinr(k,stream) = desired/max(interference+outputNoise,eps);
    end
end
end
