function resultTable = run_adaptive_mimo_mode(cfg, channelName)
%RUN_ADAPTIVE_MIMO_MODE Joint adaptation of MIMO mode and QAM order.
% The link estimates the channel from orthogonal pilots, predicts uncoded QAM
% BER for spatial-MMSE and Alamouti combining, then chooses the highest net
% spectral-efficiency candidate predicted to meet targetBER. If none qualifies,
% it selects the candidate with the lowest predicted BER and flags an outage.
%
% The channel is block-static over the complete frame and both STBC slots.
% Reported goodput is bit-level correctly received payload rate, not CRC-goodput.

if ~isfield(cfg,'adaptiveMimo')
    error('run_adaptive_mimo_mode:MissingConfig','cfg.adaptiveMimo is required.');
end
required = {'snrDb','modulations','numFrames','targetBER'};
for q = 1:numel(required)
    if ~isfield(cfg.adaptiveMimo,required{q})
        error('run_adaptive_mimo_mode:MissingField','Missing cfg.adaptiveMimo.%s.',required{q});
    end
end
if cfg.system.nTx ~= 2 || mod(cfg.system.nDataSymbols,2) ~= 0
    error('run_adaptive_mimo_mode:Topology', ...
        'Joint mode adaptation currently requires two Tx antennas and an even data-symbol count.');
end

N = cfg.system.nFFT;
Nt = cfg.system.nTx;
Nr = cfg.system.nRx;
Ns = cfg.system.nDataSymbols;
cp = cfg.system.cpLength;
Ppilot = cfg.system.pilotPower;
mods = cfg.adaptiveMimo.modulations(:).';
numFrames = cfg.adaptiveMimo.numFrames;
targetBER = cfg.adaptiveMimo.targetBER;
if isempty(mods) || any(diff(mods) <= 0) || numFrames < 1 || numFrames ~= floor(numFrames) || ...
        any(mods < 4) || any(abs(log2(mods)-round(log2(mods))) > eps) || ...
        any(abs(log2(sqrt(mods))-round(log2(sqrt(mods)))) > eps)
    error('run_adaptive_mimo_mode:InvalidSweep','Invalid modulation list or frame count.');
end
if ~isscalar(targetBER) || ~isfinite(targetBER) || targetBER <= 0 || targetBER >= 0.5
    error('run_adaptive_mimo_mode:InvalidTarget','targetBER must lie in (0,0.5).');
end
if isfield(cfg.adaptiveMimo,'bootstrapSamples')
    bootstrapSamples = cfg.adaptiveMimo.bootstrapSamples;
else
    bootstrapSamples = cfg.statistics.bootstrapSamples;
end
if isfield(cfg.adaptiveMimo,'confidenceLevel')
    confidenceLevel = cfg.adaptiveMimo.confidenceLevel;
else
    confidenceLevel = cfg.statistics.confidenceLevel;
end
if isfield(cfg.adaptiveMimo,'bootstrapSeed')
    bootstrapSeed = cfg.adaptiveMimo.bootstrapSeed;
else
    bootstrapSeed = cfg.statistics.bootstrapSeed + 9000;
end
if isfield(cfg.adaptiveMimo,'frameSeedBase') && ~isempty(cfg.adaptiveMimo.frameSeedBase)
    frameSeedBase = cfg.adaptiveMimo.frameSeedBase;
    if ~isscalar(frameSeedBase) || ~isreal(frameSeedBase) || ...
            ~isfinite(frameSeedBase) || frameSeedBase < 0 || frameSeedBase ~= floor(frameSeedBase)
        error('run_adaptive_mimo_mode:FrameSeed','frameSeedBase must be a nonnegative integer scalar.');
    end
    useFrameSeedSchedule = true;
else
    frameSeedBase = 0;
    useFrameSeedSchedule = false;
end

[meanChannel,covarianceScale] = channel_estimation_prior(cfg,channelName);
useSpatialAwareLMMSE = isfield(cfg,'estimation') && ...
    isfield(cfg.estimation,'spatialAwareLMMSE') && logical(cfg.estimation.spatialAwareLMMSE);
if useSpatialAwareLMMSE
    [~,~,RrxPrior,RtxPrior] = spatial_correlation_factors(cfg);
else
    RtxPrior = eye(Nt);
    RrxPrior = eye(Nr);
end
resultTable = [];
fprintf('\n[Adaptive MIMO mode + QAM] Channel: %s\n',upper(char(channelName)));

for snrDb = cfg.adaptiveMimo.snrDb
    noiseVarianceFD = 10^(-snrDb/10);
    noiseVarianceTime = noiseVarianceFD/N;
    if useSpatialAwareLMMSE
        jointModel = channel_estimation_joint_lmmse_setup(cfg.channel.pdp, ...
            noiseVarianceFD,Ppilot,N,covarianceScale,RtxPrior,RrxPrior);
        channelErrorVariance = jointModel.posteriorErrorVarianceByTone;
        % Retain the full receive-side covariance induced by channel-estimation
        % error. The diagonal-only scalar approximation can hide the benefit
        % of spatially-aware estimation at the detector input.
        receiverNoiseInput = channel_error_effective_noise_covariance( ...
            jointModel,noiseVarianceFD,1);
        receiverNoiseModel = "full-spatial-covariance";
    else
        channelErrorVariance = channel_estimation_error_variance( ...
            cfg.channel.pdp,noiseVarianceFD,Ppilot,N,covarianceScale);
        receiverNoiseInput = noiseVarianceFD + Nt*channelErrorVariance;
        receiverNoiseModel = "scalar-per-tone";
    end

    totalBits = 0;
    totalBitErrors = 0;
    totalFrameErrors = 0;
    totalCorrectBits = 0;
    frameBitErrors = zeros(numFrames,1);
    framePayloadBits = zeros(numFrames,1);
    framePredictedBER = zeros(numFrames,1);
    frameTargetFeasible = false(numFrames,1);
    frameModeIndex = zeros(numFrames,1);
    frameModulation = zeros(numFrames,1);
    frameSpectralEfficiency = zeros(numFrames,1);
    modeCounts = zeros(1,2);
    modulationCounts = zeros(1,numel(mods));

    for frame = 1:numFrames
        % Optional frame-level reseeding keeps channel and pilot-noise
        % realizations paired across estimators even if their previous-frame
        % mode/modulation choices consumed different numbers of random values.
        if useFrameSeedSchedule
            rng(frameSeedBase+frame,'twister');
        end
        % One channel realization and one channel estimate are used to score
        % both candidate architectures before the selected architecture sends data.
        taps = generate_mimo_channel(cfg,channelName);
        pilotGrid = complex(zeros(N,Nt,Nt));
        for tx = 1:Nt
            pilotGrid(:,tx,tx) = sqrt(Ppilot);
        end
        pilotTx = ofdm_modulate(pilotGrid,cp);
        pilotRx = apply_mimo_ofdm_channel(pilotTx,taps,noiseVarianceTime);
        pilotY = ofdm_demodulate(pilotRx,N,cp);
        Hls = channel_estimation_ls(pilotY,Ppilot);
        if useSpatialAwareLMMSE
            Hest = channel_estimation_joint_lmmse_apply(Hls,meanChannel,jointModel);
        else
            Hest = channel_estimation_lmmse(Hls,cfg.channel.pdp,noiseVarianceFD, ...
                Ppilot,meanChannel,covarianceScale);
        end

        % Candidate 1: two spatial streams with linear MMSE detection.
        [~,spatialSinr] = mimo_equalize(complex(zeros(N,Nr,1)),Hest, ...
            receiverNoiseInput,'mmse');
        % Candidate 2: rate-one Alamouti with coherent orthogonal combining.
        [~,alamoutiSinr] = alamouti_decode(complex(zeros(N,Nr,Ns)),Hest, ...
            receiverNoiseInput);
        sinrByMode = {spatialSinr(:),alamoutiSinr(:)};

        predictedBER = zeros(2,numel(mods));
        spectralEfficiency = zeros(2,numel(mods));
        for modeIndex = 1:2
            for modulationIndex = 1:numel(mods)
                predictedBER(modeIndex,modulationIndex) = mean( ...
                    qam_ber_approx(sinrByMode{modeIndex},mods(modulationIndex)));
                if modeIndex == 1
                    modeName = 'spatial-multiplexing';
                else
                    modeName = 'alamouti-stbc';
                end
                spectralEfficiency(modeIndex,modulationIndex) = ...
                    mimo_mode_spectral_efficiency(cfg,mods(modulationIndex),modeName);
            end
        end
        choice = select_mimo_link_mode(predictedBER,mods,targetBER,spectralEfficiency);

        selectedM = choice.ModulationOrder;
        bitsPerSymbol = log2(selectedM);
        if choice.ModeIndex == 1
            txBits = randi([0 1],N*Nt*Ns*bitsPerSymbol,1);
            Xdata = reshape(qam_mod(txBits,selectedM),N,Nt,Ns);
            dataTx = ofdm_modulate(Xdata,cp);
            dataRx = apply_mimo_ofdm_channel(dataTx,taps,noiseVarianceTime);
            dataY = ofdm_demodulate(dataRx,N,cp);
            [Xhat,~] = mimo_equalize(dataY,Hest,receiverNoiseInput,'mmse');
            rxBits = qam_demod(Xhat(:),selectedM);
        else
            txBits = randi([0 1],N*Ns*bitsPerSymbol,1);
            Xinformation = reshape(qam_mod(txBits,selectedM),N,2,Ns/2);
            Xstbc = alamouti_encode(Xinformation);
            dataTx = ofdm_modulate(Xstbc,cp);
            dataRx = apply_mimo_ofdm_channel(dataTx,taps,noiseVarianceTime);
            dataY = ofdm_demodulate(dataRx,N,cp);
            [Xhat,~] = alamouti_decode(dataY,Hest,receiverNoiseInput);
            rxBits = qam_demod(Xhat(:),selectedM);
        end

        errors = sum(rxBits ~= txBits);
        nBits = numel(txBits);
        totalBits = totalBits+nBits;
        totalBitErrors = totalBitErrors+errors;
        totalFrameErrors = totalFrameErrors+double(errors > 0);
        totalCorrectBits = totalCorrectBits+(nBits-errors);
        frameBitErrors(frame) = errors;
        framePayloadBits(frame) = nBits;
        framePredictedBER(frame) = choice.PredictedBER;
        frameTargetFeasible(frame) = choice.TargetFeasible;
        frameModeIndex(frame) = choice.ModeIndex;
        frameModulation(frame) = selectedM;
        frameSpectralEfficiency(frame) = choice.SpectralEfficiency_bpsHz;
        modeCounts(choice.ModeIndex) = modeCounts(choice.ModeIndex)+1;
        modulationIndex = find(mods == selectedM,1);
        modulationCounts(modulationIndex) = modulationCounts(modulationIndex)+1;
    end

    metrics = bit_error_metrics(totalBitErrors,totalBits,totalFrameErrors,numFrames);
    ferMetrics = bit_error_metrics(totalFrameErrors,numFrames,0,1);
    seed = bootstrapSeed + round(snrDb*100) + sum(double(char(channelName)));
    clusterCI = cluster_bootstrap_ber_ci(frameBitErrors,framePayloadBits, ...
        bootstrapSamples,confidenceLevel,seed);
    predictedBERFrameAverage = mean(framePredictedBER);
    predictedBERPayloadWeighted = sum(framePredictedBER.*framePayloadBits)/sum(framePayloadBits);
    predictedFeasibleRatio = mean(frameTargetFeasible);
    fallbackRatio = 1-predictedFeasibleRatio;
    observedTargetSatisfied = metrics.BER <= targetBER;
    wilsonUpperBelowTarget = metrics.BER_CI95_Upper <= targetBER;
    bootstrapUpperBelowTarget = clusterCI(2) <= targetBER;
    bothIntervalsBelowTarget = wilsonUpperBelowTarget && bootstrapUpperBelowTarget;
    spatialRatio = modeCounts(1)/numFrames;
    alamoutiRatio = modeCounts(2)/numFrames;
    spectralEfficiencyMean = mean(frameSpectralEfficiency);
    goodput = (totalCorrectBits/numFrames)/((Ns+Nt)*N)*N/(N+cp);
    if ndims(receiverNoiseInput) == 3
        offDiagonalRatio = zeros(N,1);
        for tone = 1:N
            Rtone = receiverNoiseInput(:,:,tone);
            offDiagonalRatio(tone) = norm(Rtone-diag(diag(Rtone)),'fro') / ...
                max(real(trace(Rtone)),eps);
        end
        meanOffDiagonalRatio = mean(offDiagonalRatio);
    else
        meanOffDiagonalRatio = 0;
    end

    row = table(string(channelName),string(receiverNoiseModel),snrDb,metrics.BER,metrics.BLER, ...
        metrics.BER_CI95_Lower,metrics.BER_CI95_Upper, ...
        ferMetrics.BER_CI95_Lower,ferMetrics.BER_CI95_Upper, ...
        clusterCI(1),clusterCI(2),spectralEfficiencyMean,goodput, ...
        predictedBERFrameAverage,predictedBERPayloadWeighted,predictedFeasibleRatio, ...
        fallbackRatio,spatialRatio,alamoutiRatio, ...
        modulationCounts(1)/numFrames,modulationCounts(min(2,numel(mods)))/numFrames, ...
        modulationCounts(end)/numFrames,observedTargetSatisfied,wilsonUpperBelowTarget, ...
        bootstrapUpperBelowTarget,bothIntervalsBelowTarget,totalBits,numFrames, ...
        mean(channelErrorVariance),meanOffDiagonalRatio, ...
        'VariableNames',{'Channel','ReceiverNoiseModel','EsN0_dB','BER','FER','BER_Wilson_CI95_Lower', ...
        'BER_Wilson_CI95_Upper','FER_Wilson_CI95_Lower','FER_Wilson_CI95_Upper', ...
        'BER_FrameBootstrap_CI95_Lower','BER_FrameBootstrap_CI95_Upper', ...
        'SpectralEfficiency_bpsHz','BitGoodput_bpsHz','MeanPredictedBER_FrameAverage', ...
        'PredictedBER_PayloadWeighted','PredictedTargetFeasibleFrameRatio', ...
        'FallbackFrameRatio','SpatialMultiplexingFrameRatio','AlamoutiFrameRatio', ...
        'QPSK_FrameRatio','16QAM_FrameRatio','64QAM_FrameRatio', ...
        'ObservedBER_TargetSatisfied','ObservedBER_WilsonCI95UpperBelowTarget', ...
        'ObservedBER_ClusterBootstrapCI95UpperBelowTarget', ...
        'ObservedBER_BothCI95UpperBelowTarget','PayloadBits','NumFrames', ...
        'MeanChannelEstimationErrorVariance','MeanEffectiveNoiseOffDiagonalRatio'});
    if isempty(resultTable)
        resultTable = row;
    else
        resultTable = [resultTable;row]; %#ok<AGROW>
    end

    fprintf(['  Es/N0=%2d dB | BER=%.3g [Wilson %.2g, %.2g; bootstrap %.2g, %.2g] | ' ...
        'SE=%.2f | mode mix=[SM %.0f%%, Alamouti %.0f%%] | feasible=%.0f%%\n'], ...
        snrDb,metrics.BER,metrics.BER_CI95_Lower,metrics.BER_CI95_Upper, ...
        clusterCI(1),clusterCI(2),spectralEfficiencyMean,100*spatialRatio, ...
        100*alamoutiRatio,100*predictedFeasibleRatio);
end
end
