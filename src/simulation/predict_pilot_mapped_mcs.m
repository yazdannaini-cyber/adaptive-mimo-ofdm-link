function prediction = predict_pilot_mapped_mcs(Hestimate, noiseVariance, modulations, targetBER, spectralEfficiency)
%PREDICT_PILOT_MAPPED_MCS Score SM and Alamouti/QAM candidates from pilot CSI.
% The channel estimate is the causal prediction at the first data OFDM symbol.
% Candidate BER uses the project's uncoded Gray-QAM AWGN approximation applied
% to post-equalization SINR. It is a link-adaptation predictor, not a guarantee.

if ~isnumeric(Hestimate) || ndims(Hestimate) > 3 || size(Hestimate,3) ~= 2 || ...
        isempty(Hestimate) || any(~isfinite(real(Hestimate(:)))) || any(~isfinite(imag(Hestimate(:))))
    error('predict_pilot_mapped_mcs:InvalidChannel', ...
        'Hestimate must be finite [nFFT x nRx x 2] data.');
end
if ~isscalar(noiseVariance) || ~isreal(noiseVariance) || ...
        ~isfinite(noiseVariance) || noiseVariance <= 0
    error('predict_pilot_mapped_mcs:InvalidNoise', ...
        'noiseVariance must be a finite positive scalar.');
end
modulations = double(modulations(:).');
if isempty(modulations) || any(diff(modulations) <= 0) || ...
        any(modulations < 4) || any(abs(log2(modulations)-round(log2(modulations))) > eps) || ...
        any(abs(log2(sqrt(modulations))-round(log2(sqrt(modulations)))) > eps)
    error('predict_pilot_mapped_mcs:InvalidModulations', ...
        'modulations must be ascending square power-of-four QAM orders.');
end
if ~isscalar(targetBER) || ~isreal(targetBER) || ~isfinite(targetBER) || ...
        targetBER <= 0 || targetBER >= 0.5
    error('predict_pilot_mapped_mcs:InvalidTarget', 'targetBER must lie in (0,0.5).');
end
if ~isnumeric(spectralEfficiency) || ~isequal(size(spectralEfficiency),[2 numel(modulations)]) || ...
        any(~isfinite(spectralEfficiency(:))) || any(spectralEfficiency(:) <= 0)
    error('predict_pilot_mapped_mcs:InvalidEfficiency', ...
        'spectralEfficiency must be a positive 2-by-numel(modulations) matrix.');
end

[nFFT,nRx,~] = size(Hestimate);
[~,smPostSinr] = mimo_equalize(complex(zeros(nFFT,nRx,1)),Hestimate,noiseVariance,'mmse');
alamoutiPostSinr = alamouti_timevarying_postsinr(Hestimate,Hestimate,noiseVariance);
predictedBER = zeros(2,numel(modulations));
for m = 1:numel(modulations)
    predictedBER(1,m) = mean(qam_ber_approx(smPostSinr(:),modulations(m)));
    predictedBER(2,m) = mean(qam_ber_approx(alamoutiPostSinr(:),modulations(m)));
end
choice = select_mimo_link_mode(predictedBER,modulations,targetBER,spectralEfficiency);

prediction.predictedBER = predictedBER;
prediction.spectralEfficiency = spectralEfficiency;
prediction.modulations = modulations;
prediction.postSinrSpatialMultiplexing = smPostSinr;
prediction.postSinrAlamouti = alamoutiPostSinr;
prediction.choice = choice;
end

function postSinr = alamouti_timevarying_postsinr(Hfirst,Hsecond,noiseVariance)
% Candidate-side version of the generalized Alamouti effective channel.
[nFFT,nRx,~] = size(Hfirst);
postSinr = zeros(nFFT,2);
regularization = max(double(noiseVariance),1e-12);
for k = 1:nFFT
    H1 = reshape(Hfirst(k,:,:),nRx,2);
    H2 = reshape(Hsecond(k,:,:),nRx,2);
    A = [H1(:,1), H1(:,2); conj(H2(:,2)), -conj(H2(:,1))];
    W = (A'*A + regularization*eye(2))\A';
    G = W*A;
    for stream = 1:2
        desired = abs(G(stream,stream))^2;
        interference = max(sum(abs(G(stream,:)).^2)-desired,0);
        outputNoise = double(noiseVariance)*sum(abs(W(stream,:)).^2);
        postSinr(k,stream) = desired/max(interference+outputNoise,eps);
    end
end
end
