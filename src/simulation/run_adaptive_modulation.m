function resultTable = run_adaptive_modulation(cfg, channelName)
%RUN_ADAPTIVE_MODULATION Select QAM order using predicted post-SINR BER.
% SNR axis is Es/N0 per subcarrier, not Eb/N0, because QAM order changes.
% The predictor includes posterior channel-estimation uncertainty.

fprintf('\n[Adaptive modulation] Channel: %s\n',upper(channelName));
resultTable = [];
N = cfg.system.nFFT;
Nt = cfg.system.nTx;
Nr = cfg.system.nRx;
Ns = cfg.system.nDataSymbols;
cp = cfg.system.cpLength;
Ppilot = cfg.system.pilotPower;
numFrames = cfg.simulation.numFrames;
mods = cfg.adaptive.modulations;
thresholds = cfg.adaptive.thresholdDb;
selectionMethod = lower(char(cfg.adaptive.selectionMethod));
targetBER = cfg.adaptive.targetBER;
[meanChannel,covarianceScale] = channel_estimation_prior(cfg,channelName);

if numel(mods) ~= numel(thresholds) || any(diff(thresholds) <= 0)
    error('run_adaptive_modulation:ThresholdConfig','Each modulation needs an increasing SINR threshold.');
end
if ~strcmp(selectionMethod,'predicted-ber') && ~strcmp(selectionMethod,'sinr-threshold')
    error('run_adaptive_modulation:SelectionMethod', ...
        'selectionMethod must be ''predicted-ber'' or ''sinr-threshold''.');
end

for snrDb = cfg.adaptive.snrDb
    noiseVarianceFD = 10^(-snrDb/10);
    noiseVarianceTime = noiseVarianceFD/N;
    totalBits = 0;
    bitErrors = 0;
    blockErrors = 0;
    totalCorrectBits = 0;
    selectedCounts = zeros(1,numel(mods));
    sinrMetricSum = 0;
    posteriorVarianceSum = 0;
    framePredictedBER = zeros(numFrames,1);
    frameBitErrors = zeros(numFrames,1);
    framePayloadBits = zeros(numFrames,1);
    frameTargetFeasible = false(numFrames,1);
    frameFallback = false(numFrames,1);

    % Posterior error variance of each estimated MIMO channel coefficient.
    % Independent transmit-stream errors add as effective noise at each Rx.
    channelErrorVariance = channel_estimation_error_variance( ...
        cfg.channel.pdp,noiseVarianceFD,Ppilot,N,covarianceScale);
    effectiveNoiseVariance = noiseVarianceFD + Nt*channelErrorVariance;
    posteriorVarianceSum = sum(channelErrorVariance);

    for frame = 1:numFrames
        taps = generate_mimo_channel(cfg,channelName);
        Xpilot = complex(zeros(N,Nt,Nt));
        for tx = 1:Nt
            Xpilot(:,tx,tx) = sqrt(Ppilot);
        end
        pilotTx = ofdm_modulate(Xpilot,cp);
        pilotRx = apply_mimo_ofdm_channel(pilotTx,taps,noiseVarianceTime);
        pilotYfd = ofdm_demodulate(pilotRx,N,cp);
        Hls = channel_estimation_ls(pilotYfd,Ppilot);
        Hlmmse = channel_estimation_lmmse(Hls,cfg.channel.pdp,noiseVarianceFD, ...
            Ppilot,meanChannel,covarianceScale);

        [~,sinrValues] = mimo_equalize(zeros(N,Nr,1),Hlmmse, ...
            effectiveNoiseVariance,'mmse');
        sortedSinr = sort(sinrValues(:));
        percentileIndex = max(1,ceil(cfg.adaptive.percentile/100*numel(sortedSinr)));
        robustSinrDb = 10*log10(max(sortedSinr(percentileIndex),eps));
        sinrMetricSum = sinrMetricSum+robustSinrDb;

        predictedBERByMode = zeros(1,numel(mods));
        for modeIndex = 1:numel(mods)
            predictedBERByMode(modeIndex) = mean(qam_ber_approx(sinrValues(:),mods(modeIndex)));
        end
        if strcmp(selectionMethod,'predicted-ber')
            [selectedIndex, targetFeasible] = ...
                select_qam_by_predicted_ber(predictedBERByMode,targetBER);
            fallback = ~targetFeasible;
        else
            eligible = find(robustSinrDb >= thresholds,1,'last');
            fallback = isempty(eligible);
            if fallback
                selectedIndex = 1; % QPSK minimum-service fallback
            else
                selectedIndex = eligible;
            end
            targetFeasible = predictedBERByMode(selectedIndex) <= targetBER;
        end

        M = mods(selectedIndex);
        selectedPredictedBER = predictedBERByMode(selectedIndex);
        bitsPerSymbol = log2(M);
        txBits = randi([0 1],N*Nt*Ns*bitsPerSymbol,1);
        Xdata = reshape(qam_mod(txBits,M),N,Nt,Ns);
        dataTx = ofdm_modulate(Xdata,cp);
        dataRx = apply_mimo_ofdm_channel(dataTx,taps,noiseVarianceTime);
        Ydata = ofdm_demodulate(dataRx,N,cp);
        % Include channel-estimation uncertainty in MMSE regularization.
        [Xhat,~] = mimo_equalize(Ydata,Hlmmse,effectiveNoiseVariance,'mmse');
        rxBits = qam_demod(Xhat(:),M);
        errors = sum(rxBits ~= txBits);
        nBits = numel(txBits);
        totalBits = totalBits+nBits;
        bitErrors = bitErrors+errors;
        blockErrors = blockErrors+double(errors > 0);
        totalCorrectBits = totalCorrectBits+(nBits-errors);

        framePredictedBER(frame) = selectedPredictedBER;
        frameBitErrors(frame) = errors;
        framePayloadBits(frame) = nBits;
        frameTargetFeasible(frame) = targetFeasible;
        frameFallback(frame) = fallback;
        selectedCounts(selectedIndex) = selectedCounts(selectedIndex)+1;
    end

    m = bit_error_metrics(bitErrors,totalBits,blockErrors,numFrames);
    adaptiveStats = adaptive_link_metrics(framePredictedBER,framePayloadBits, ...
        frameTargetFeasible,frameFallback,m.BER,m.BER_CI95_Upper,targetBER);
    clusterCI = cluster_bootstrap_ber_ci(frameBitErrors,framePayloadBits, ...
        cfg.statistics.bootstrapSamples,cfg.statistics.confidenceLevel, ...
        cfg.statistics.bootstrapSeed + round(snrDb*100) + sum(double(channelName)));
    averageK = sum(selectedCounts.*log2(mods))/numFrames;
    spectralEfficiency = Nt*averageK*Ns/(Ns+Nt)*N/(N+cp);
    goodput = (totalCorrectBits/numFrames)/((Ns+Nt)*N)*N/(N+cp);
    meanPosteriorVariance = posteriorVarianceSum/N;
    row = table(string(channelName),snrDb,m.BER,m.BLER,m.BER_CI95_Lower,m.BER_CI95_Upper, ...
        spectralEfficiency,goodput,sinrMetricSum/numFrames, ...
        selectedCounts(1)/numFrames,selectedCounts(min(2,numel(mods)))/numFrames, ...
        selectedCounts(end)/numFrames,totalBits,numFrames,targetBER, ...
        adaptiveStats.MeanPredictedBER_FrameAverage, ...
        adaptiveStats.PredictedBER_PayloadWeighted, ...
        adaptiveStats.PredictedTargetMetFrameRatio, ...
        adaptiveStats.FallbackQPSKFrameRatio, ...
        adaptiveStats.ObservedBER_TargetSatisfied, ...
        adaptiveStats.ObservedBER_CI95UpperBelowTarget, ...
        clusterCI(1) <= targetBER, ...
        (adaptiveStats.ObservedBER_CI95UpperBelowTarget && clusterCI(2) <= targetBER), ...
        clusterCI(1),clusterCI(2),meanPosteriorVariance, ...
        'VariableNames',{'Channel','EsN0_dB','BER','BLER','BER_CI95_Lower','BER_CI95_Upper', ...
        'SpectralEfficiency_bpsHz','Goodput_bpsHz','MeanRobustSINR_dB', ...
        'QPSK_FrameRatio','16QAM_FrameRatio','64QAM_FrameRatio','PayloadBits','NumFrames', ...
        'TargetBER','MeanPredictedBER_FrameAverage','PredictedBER_PayloadWeighted', ...
        'PredictedTargetMetFrameRatio','FallbackQPSKFrameRatio', ...
        'ObservedBER_TargetSatisfied','ObservedBER_WilsonCI95UpperBelowTarget', ...
        'ObservedBER_ClusterBootstrapCI95UpperBelowTarget', ...
        'ObservedBER_BothCI95UpperBelowTarget', ...
        'BER_ClusterBootstrap_CI95_Lower','BER_ClusterBootstrap_CI95_Upper', ...
        'MeanChannelEstimationErrorVariance'});
    if isempty(resultTable)
        resultTable = row;
    else
        resultTable = [resultTable;row]; %#ok<AGROW>
    end
    fprintf(['  Es/N0=%2d dB | BER=%.3g [Wilson %.2g, %.2g; frame-bootstrap %.2g, %.2g] | predBER(bit-weighted)=%.3g | ' ...
        'SE=%.2f bit/s/Hz | predicted-feasible=%.0f%% | QPSK-fallback=%.0f%% | modes=[%s]\n'], ...
        snrDb,m.BER,m.BER_CI95_Lower,m.BER_CI95_Upper,clusterCI(1),clusterCI(2), ...
        adaptiveStats.PredictedBER_PayloadWeighted,spectralEfficiency, ...
        100*adaptiveStats.PredictedTargetMetFrameRatio, ...
        100*adaptiveStats.FallbackQPSKFrameRatio,sprintf(' %.0f%%',100*selectedCounts/numFrames));
end
end
