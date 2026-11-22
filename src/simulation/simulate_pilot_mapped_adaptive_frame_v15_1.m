function frame = simulate_pilot_mapped_adaptive_frame_v15_1(cfg,modelName,params)
%SIMULATE_PILOT_MAPPED_ADAPTIVE_FRAME_V15_1 Block-adaptive pilot-mapped MIMO-OFDM.
% V15.1 improvements:
%  1) A fresh MIMO-mode/QAM decision is made after each orthogonal pilot occasion.
%  2) The decision uses only pilots observed by that occasion (strictly causal).
%  3) Per-symbol temporal-LMMSE posterior variance is added to effective noise.
%  4) Residual ICI is estimated from causal pilot innovations, with a small-
%     Doppler analytical floor at the first pilot occasion.
%  5) SM and generalized two-slot Alamouti predictions use the predicted
%     channel at the actual data symbols in the upcoming burst.
%
% References remain fixed SM-16QAM and fixed Alamouti-QPSK. All waveforms use
% the same sample-varying channel, pilot grid and receiver thermal-noise array.

required = {'nSymbols','pilotSpacing','pilotSNRdB','maxDopplerHz','txRho','rxRho', ...
    'dataSeed','pilotSeed','channelSeed','noiseSeed','targetBER','modulations'};
if ~isstruct(params) || ~all(isfield(params,required))
    error('simulate_pilot_mapped_adaptive_frame_v15_1:InvalidParams', ...
        'params is missing one or more required fields.');
end
if ~isstruct(cfg) || ~isfield(cfg,'system') || ~isfield(cfg,'channel')
    error('simulate_pilot_mapped_adaptive_frame_v15_1:InvalidConfig','Invalid project configuration.');
end
if cfg.system.nTx ~= 2 || cfg.system.nRx < 2
    error('simulate_pilot_mapped_adaptive_frame_v15_1:Topology', ...
        'V15.1 requires two transmit antennas and at least two receive antennas.');
end
if ~isscalar(params.nSymbols) || ~isnumeric(params.nSymbols) || ~isreal(params.nSymbols) || ...
        ~isfinite(params.nSymbols) || params.nSymbols ~= floor(params.nSymbols) || ...
        ~isscalar(params.pilotSpacing) || ~isnumeric(params.pilotSpacing) || ...
        ~isreal(params.pilotSpacing) || ~isfinite(params.pilotSpacing) || ...
        params.pilotSpacing ~= floor(params.pilotSpacing) || ...
        params.nSymbols < 4 || params.pilotSpacing < cfg.system.nTx
    error('simulate_pilot_mapped_adaptive_frame_v15_1:Schedule', ...
        'nSymbols and pilotSpacing must be valid integers.');
end
if ~isscalar(params.pilotSNRdB) || ~isreal(params.pilotSNRdB) || ~isfinite(params.pilotSNRdB) || ...
        ~isscalar(params.maxDopplerHz) || ~isreal(params.maxDopplerHz) || ...
        ~isfinite(params.maxDopplerHz) || params.maxDopplerHz < 0
    error('simulate_pilot_mapped_adaptive_frame_v15_1:LinkParameters', ...
        'Pilot SNR must be finite and Doppler must be nonnegative.');
end
if ~isscalar(params.txRho) || ~isscalar(params.rxRho) || ...
        ~isreal(params.txRho) || ~isreal(params.rxRho) || ...
        ~isfinite(params.txRho) || ~isfinite(params.rxRho) || ...
        params.txRho < 0 || params.rxRho < 0 || params.txRho > 0.99 || params.rxRho > 0.99
    error('simulate_pilot_mapped_adaptive_frame_v15_1:Correlation', ...
        'Spatial-correlation coefficients must lie in [0,0.99].');
end
if ~isscalar(params.targetBER) || ~isreal(params.targetBER) || ...
        ~isfinite(params.targetBER) || params.targetBER <= 0 || params.targetBER >= 0.5
    error('simulate_pilot_mapped_adaptive_frame_v15_1:TargetBER', ...
        'targetBER must lie in (0,0.5).');
end
modulations = double(params.modulations(:).');
if isempty(modulations) || ~all(ismember([4 16],modulations)) || ...
        any(diff(modulations) <= 0) || ...
        any(abs(log2(modulations)-round(log2(modulations))) > eps) || ...
        any(abs(log2(sqrt(modulations))-round(log2(sqrt(modulations)))) > eps)
    error('simulate_pilot_mapped_adaptive_frame_v15_1:Modulations', ...
        'modulations must be ascending square power-of-four QAM orders including 4 and 16.');
end
seeds = [params.dataSeed params.pilotSeed params.channelSeed params.noiseSeed];
if any(~isfinite(seeds)) || any(seeds < 0) || any(seeds ~= floor(seeds))
    error('simulate_pilot_mapped_adaptive_frame_v15_1:Seeds','All seeds must be nonnegative integers.');
end
if isfield(params,'estimator') && ~strcmpi(char(params.estimator),'jakes-lmmse')
    error('simulate_pilot_mapped_adaptive_frame_v15_1:Estimator', ...
        'V15.1 currently uses the causal Jakes-LMMSE tracker.');
end

N = cfg.system.nFFT; cp = cfg.system.cpLength;
Nt = cfg.system.nTx; Nr = cfg.system.nRx; ns = params.nSymbols;
schedule = mimo_ofdm_pilot_schedule(ns,Nt,params.pilotSpacing);
blockDataIndices = build_data_blocks(schedule,Nt,ns);

cfgFrame = cfg;
cfgFrame.channel.doppler.maxDopplerHz = params.maxDopplerHz;
cfgFrame.channel.spatialCorrelation.enabled = true;
cfgFrame.channel.spatialCorrelation.txRho = params.txRho;
cfgFrame.channel.spatialCorrelation.rxRho = params.rxRho;

% Build the orthogonal full-band TDM pilot waveform.
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(params.pilotSeed,'twister');
nOccasions = numel(schedule.pilotOccasionStarts);
pilotSymbols = complex(zeros(N,Nt,nOccasions));
Xpilot = complex(zeros(N,Nt,ns));
for tx = 1:Nt
    for q = 1:nOccasions
        bits = randi([0 1],2*N,1);
        knownSequence = sqrt(cfg.system.pilotPower)*qam_mod(bits,4);
        pilotSymbols(:,tx,q) = knownSequence;
        Xpilot(:,tx,schedule.pilotIndicesByTx{tx}(q)) = knownSequence;
    end
end

blockLength = N+cp;
sampleRateHz = N*cfg.system.subcarrierSpacingHz;
totalSamples = blockLength*ns;
[tapSequence,components,channelMeta] = generate_time_varying_mimo_channel( ...
    cfgFrame,modelName,totalSamples,sampleRateHz,params.channelSeed);
[diffuseVariance,meanChannelPower] = channel_variance_config(cfg,modelName);
noiseVarianceFD = meanChannelPower*10^(-params.pilotSNRdB/10);
rng(params.noiseSeed,'twister');
noiseFD = sqrt(noiseVarianceFD/2)*(randn(N,Nr,ns)+1i*randn(N,Nr,ns));

% Pilot-only observation is sufficient because CP removes inter-symbol
% convolution; it still contains within-symbol Doppler-induced ICI.
pilotRx = apply_time_varying_mimo_ofdm_channel(ofdm_modulate(Xpilot,cp),tapSequence,0);
pilotY = ofdm_demodulate(pilotRx,N,cp)+noiseFD;
observationsByTx = extract_mimo_pilot_observations(pilotY,pilotSymbols,schedule);
Hmean = time_varying_ofdm_diagonal_response(components.los,N,cp,ns);
Tsym = blockLength/(N*cfg.system.subcarrierSpacingHz);

% Form causal estimates for every upcoming data burst. The pilot-noise
% covariance used by the tracker is inflated by the ICI innovation estimate
% from the current pilot occasion; the current estimate is based only on
% earlier pilots plus the current pilot innovation.
blocks = repmat(struct(),numel(blockDataIndices),1);
iciNoiseByOccasion = zeros(1,numel(blockDataIndices));
for b = 1:numel(blockDataIndices)
    dataIdx = blockDataIndices{b};
    previousHistory = iciNoiseByOccasion(1:b-1);
    if isempty(previousHistory), provisionalHistory = 0; else, provisionalHistory = [previousHistory previousHistory(end)]; end
    [Hprovisional,~] = track_available_pilots_v15_1( ...
        observationsByTx,schedule,b,Hmean,provisionalHistory, ...
        cfg,params,ns,diffuseVariance,Tsym,dataIdx,noiseVarianceFD);
    [iciVariance,iciDetail] = estimate_pilot_ici_variance( ...
        observationsByTx,schedule,b,Hmean,noiseVarianceFD, ...
        cfg.system.pilotPower,params.maxDopplerHz,Tsym,diffuseVariance, ...
        cfg.system.subcarrierSpacingHz,Hprovisional,previousHistory);
    iciNoiseByOccasion(b) = iciVariance;
    [Hburst,posteriorByTx] = track_available_pilots_v15_1( ...
        observationsByTx,schedule,b,Hmean,iciNoiseByOccasion(1:b), ...
        cfg,params,ns,diffuseVariance,Tsym,dataIdx,noiseVarianceFD);
    iciBySymbol = repmat(iciVariance,1,numel(dataIdx));
    blockSE = [Nt;1]*log2(modulations)*(numel(dataIdx)/ns)*N/(N+cp);
    prediction = predict_pilot_mapped_mcs_block(Hburst,noiseVarianceFD, ...
        modulations,params.targetBER,blockSE,posteriorByTx,iciBySymbol);

    blocks(b).dataIndices = dataIdx;
    blocks(b).HestimateForDecision = Hburst;
    blocks(b).posteriorVarianceByTx = posteriorByTx;
    blocks(b).iciNoiseVariance = iciVariance;
    blocks(b).iciDetail = iciDetail;
    blocks(b).prediction = prediction;
    blocks(b).choice = prediction.choice;
end

% Frame-level ICI variance map used by the actual receiver and fixed baselines.
iciNoiseBySymbol = zeros(1,ns);
for b = 1:numel(blocks)
    iciNoiseBySymbol(blocks(b).dataIndices) = blocks(b).iciNoiseVariance;
end

% Causal per-burst adaptive waveform and two fixed reference waveforms.
adaptiveResult = simulate_block_adaptive_candidate_v15_1(cfgFrame,params, ...
    Xpilot,pilotSymbols,schedule,blockDataIndices,blocks,tapSequence,noiseFD, ...
    noiseVarianceFD,iciNoiseBySymbol);

sm16Index = find(modulations==16,1);
qpskIndex = find(modulations==4,1);
sm16Pred = mean(arrayfun(@(x)x.prediction.predictedBER(1,sm16Index),blocks));
alQpskPred = mean(arrayfun(@(x)x.prediction.predictedBER(2,qpskIndex),blocks));
fixedSM16 = simulate_one_fixed_candidate_v15_1(cfgFrame,params, ...
    Xpilot,pilotSymbols,schedule,blockDataIndices,blocks,tapSequence,noiseFD, ...
    noiseVarianceFD,iciNoiseBySymbol,'spatial-multiplexing',16,sm16Pred,sm16Pred<=params.targetBER,1);
fixedAlamoutiQPSK = simulate_one_fixed_candidate_v15_1(cfgFrame,params, ...
    Xpilot,pilotSymbols,schedule,blockDataIndices,blocks,tapSequence,noiseFD, ...
    noiseVarianceFD,iciNoiseBySymbol,'alamouti-stbc',4,alQpskPred,alQpskPred<=params.targetBER,2);

frame.schedule = schedule;
frame.dataBlocks = blockDataIndices;
frame.channelMeta = channelMeta;
frame.noiseVarianceFD = noiseVarianceFD;
frame.targetBER = params.targetBER;
frame.estimator = 'jakes-lmmse';
frame.blocks = blocks;
frame.choice = blocks(1).choice; % backward-compatible first-burst choice
frame.prediction = blocks(1).prediction;
frame.adaptive = adaptiveResult;
frame.fixedSM16QAM = fixedSM16;
frame.fixedAlamoutiQPSK = fixedAlamoutiQPSK;
end

function [Hburst,posteriorByTx] = track_available_pilots_v15_1( ...
    observationsByTx,schedule,blockIndex,Hmean,iciNoiseByOccasion, ...
    cfg,params,nSymbols,diffuseVariance,symbolPeriodSeconds,dataIdx,thermalNoiseVariance)
N=cfg.system.nFFT; Nr=cfg.system.nRx; Nt=cfg.system.nTx;
Hburst=complex(zeros(N,Nr,Nt,numel(dataIdx)));
posteriorByTx=zeros(Nt,numel(dataIdx));
% The caller stores only the current burst predictions; use explicit slice
% assignment so MATLAB preserves the four-dimensional channel layout.
for tx=1:Nt
    pilotIndices=schedule.pilotIndicesByTx{tx}(1:blockIndex);
    obs=reshape(observationsByTx{tx}(:,:,1:blockIndex),N,Nr,1,blockIndex);
    prior=reshape(Hmean(:,:,tx,:),N,Nr,1,nSymbols);
    pilotNoiseByOccasion=pilot_ls_noise_profile(thermalNoiseVariance, ...
        cfg.system.pilotPower,iciNoiseByOccasion(1:blockIndex));
    setup=temporal_channel_tracker_setup(pilotIndices,nSymbols, ...
        pilotNoiseByOccasion,params.maxDopplerHz, ...
        symbolPeriodSeconds,diffuseVariance);
    [track,post]=temporal_channel_tracking_apply(obs,pilotIndices,nSymbols, ...
        'jakes-lmmse',prior,setup);
    Hburst(:,:,tx,:)=reshape(track(:,:,1,dataIdx),N,Nr,1,numel(dataIdx));
    posteriorByTx(tx,:)=post(dataIdx).';
end
end

function blocks = build_data_blocks(schedule,nTx,nSymbols)
starts = schedule.pilotOccasionStarts;
blocks = cell(numel(starts),1);
for b = 1:numel(starts)
    firstData = starts(b)+nTx;
    if b < numel(starts)
        lastData = starts(b+1)-1;
    else
        lastData = nSymbols;
    end
    idx = firstData:lastData;
    if isempty(idx) || mod(numel(idx),2) ~= 0
        error('simulate_pilot_mapped_adaptive_frame_v15_1:OddDataSlots', ...
            'Each pilot-to-pilot data burst must contain an even number of OFDM symbols.');
    end
    blocks{b} = idx;
end
end

function result = simulate_block_adaptive_candidate_v15_1(cfg,params, ...
    Xpilot,pilotSymbols,schedule,blockDataIndices,blocks,tapSequence,noiseFD, ...
    noiseVarianceFD,iciNoiseBySymbol)
N = cfg.system.nFFT; cp = cfg.system.cpLength;
Nt = cfg.system.nTx; Nr = cfg.system.nRx; ns = params.nSymbols;
Xfd = Xpilot;
txBitsByBlock = cell(numel(blocks),1);
modeNames = cell(numel(blocks),1);
modOrders = zeros(numel(blocks),1);
blockPredictedBER = zeros(numel(blocks),1);
blockFeasible = false(numel(blocks),1);
blockSE = zeros(numel(blocks),1);
rng(params.dataSeed+100,'twister');

for b = 1:numel(blocks)
    idx = blockDataIndices{b};
    choice = blocks(b).choice;
    modeNames{b} = choice.Mode;
    modOrders(b) = choice.ModulationOrder;
    blockPredictedBER(b) = choice.PredictedBER;
    blockFeasible(b) = choice.TargetFeasible;
    bitsPerQAM = log2(modOrders(b));
    if strcmpi(modeNames{b},'spatial-multiplexing')
        bitsBlock = zeros(N*bitsPerQAM,Nt,numel(idx));
        for q = 1:numel(idx)
            for tx = 1:Nt
                bits = randi([0 1],N*bitsPerQAM,1);
                bitsBlock(:,tx,q) = bits;
                Xfd(:,tx,idx(q)) = qam_mod(bits,modOrders(b));
            end
        end
        txBitsByBlock{b} = bitsBlock;
        blockSE(b) = Nt*bitsPerQAM;
    else
        numPairs = numel(idx)/2;
        bitsBlock = zeros(N*bitsPerQAM,2,numPairs);
        informationGrid = complex(zeros(N,2,numPairs));
        for pair = 1:numPairs
            for slot = 1:2
                bits = randi([0 1],N*bitsPerQAM,1);
                bitsBlock(:,slot,pair) = bits;
                informationGrid(:,slot,pair) = qam_mod(bits,modOrders(b));
            end
        end
        encoded = alamouti_encode(informationGrid);
        for pair = 1:numPairs
            Xfd(:,:,idx(2*pair-1)) = encoded(:,:,2*pair-1);
            Xfd(:,:,idx(2*pair)) = encoded(:,:,2*pair);
        end
        txBitsByBlock{b} = bitsBlock;
        blockSE(b) = bitsPerQAM;
    end
end

rxWaveform = apply_time_varying_mimo_ofdm_channel(ofdm_modulate(Xfd,cp),tapSequence,0);
YfdNoiseless = ofdm_demodulate(rxWaveform,N,cp);
Yfd = YfdNoiseless+noiseFD;
observations = extract_mimo_pilot_observations(Yfd,pilotSymbols,schedule);
[Hest,posteriorVariance] = assemble_burst_channel_estimates_v15_1( ...
    blocks,blockDataIndices,N,Nr,Nt,ns);

bitErrors = 0; totalBits = 0;
for b = 1:numel(blocks)
    idx = blockDataIndices{b};
    modeName = modeNames{b};
    M = modOrders(b);
    txBits = txBitsByBlock{b};
    rxBits = zeros(size(txBits));
    if strcmpi(modeName,'spatial-multiplexing')
        for q = 1:numel(idx)
            s = idx(q);
            Hs = reshape(Hest(:,:,:,s),N,Nr,Nt);
            sigmaEff = noiseVarianceFD + sum(posteriorVariance(:,s)) + iciNoiseBySymbol(s);
            [Xhat,~] = mimo_equalize(reshape(Yfd(:,:,s),N,Nr,1),Hs,sigmaEff,'mmse');
            for tx = 1:Nt
                rxBits(:,tx,q) = qam_demod(Xhat(:,tx,1),M);
            end
        end
    else
        for pair = 1:numel(idx)/2
            s1 = idx(2*pair-1); s2 = idx(2*pair);
            Ypair = cat(3,Yfd(:,:,s1),Yfd(:,:,s2));
            H1 = reshape(Hest(:,:,:,s1),N,Nr,Nt);
            H2 = reshape(Hest(:,:,:,s2),N,Nr,Nt);
            sigma1 = noiseVarianceFD + sum(posteriorVariance(:,s1)) + iciNoiseBySymbol(s1);
            sigma2 = noiseVarianceFD + sum(posteriorVariance(:,s2)) + iciNoiseBySymbol(s2);
            [Xhat,~] = alamouti_decode_timevarying(Ypair,H1,H2,[sigma1 sigma2]);
            rxBits(:,1,pair) = qam_demod(Xhat(:,1),M);
            rxBits(:,2,pair) = qam_demod(Xhat(:,2),M);
        end
    end
    bitErrors = bitErrors+nnz(rxBits(:)~=txBits(:));
    totalBits = totalBits+numel(txBits);
end

[commonMetrics,iciRatio] = measure_candidate_channel_quality_v15_1( ...
    Xfd,YfdNoiseless,Hest,posteriorVariance,tapSequence,N,cp,ns,schedule,observations);
netSE = sum(blockSE.*cellfun(@numel,blockDataIndices(:)))/ns*N/(N+cp);
result.Mode = 'block-adaptive';
if numel(unique(modeNames)) == 1 && numel(unique(modOrders)) == 1
    result.Mode = modeNames{1};
end
result.ModulationOrder = NaN;
if numel(unique(modOrders)) == 1, result.ModulationOrder = modOrders(1); end
result.ModeByBlock = modeNames;
result.ModulationOrderByBlock = modOrders;
result.SpatialMultiplexingFrameRatio = sum(strcmp(modeNames,'spatial-multiplexing').*cellfun(@numel,blockDataIndices(:)))/schedule.numDataSymbols;
result.AlamoutiFrameRatio = 1-result.SpatialMultiplexingFrameRatio;
result.QPSKFrameRatio = sum((modOrders==4).*cellfun(@numel,blockDataIndices(:)))/schedule.numDataSymbols;
result.QAM16FrameRatio = sum((modOrders==16).*cellfun(@numel,blockDataIndices(:)))/schedule.numDataSymbols;
result.PredictedBER = sum(blockPredictedBER.*cellfun(@numel,blockDataIndices(:)))/sum(cellfun(@numel,blockDataIndices));
result.TargetFeasible = all(blockFeasible);
result.TargetFeasibleBlockRatio = sum(blockFeasible.*cellfun(@numel,blockDataIndices(:)))/schedule.numDataSymbols;
result.BER = bitErrors/max(totalBits,1);
result.FER = double(bitErrors>0);
result.BitErrors = bitErrors;
result.TotalBits = totalBits;
result.NMSE = commonMetrics.NMSE;
result.MeanPosteriorVariance = commonMetrics.MeanPosteriorVariance;
result.MeanEffectiveNoiseVariance = mean(noiseVarianceFD + sum(posteriorVariance(:,schedule.dataIndices),1) + iciNoiseBySymbol(schedule.dataIndices));
result.MeanEstimatedICINoiseVariance = mean(iciNoiseBySymbol(schedule.dataIndices));
result.ICI_Ratio = iciRatio;
result.PilotContaminationVariance = commonMetrics.PilotContaminationVariance;
result.PilotOverhead = schedule.pilotOverhead;
result.NumPilotSymbols = schedule.numPilotSymbols;
result.NumDataSymbols = schedule.numDataSymbols;
result.NetSpectralEfficiency_bpsHz = netSE;
result.UncodedGoodput_bpsHz = netSE*(1-result.BER);
result.FrameError = bitErrors>0;
result.NumAdaptationBlocks = numel(blocks);
result.CausalReceiverEstimates = true;
result.BlockPredictedBER = blockPredictedBER;
result.BlockTargetFeasible = blockFeasible;
end

function result = simulate_one_fixed_candidate_v15_1(cfg,params, ...
    Xpilot,pilotSymbols,schedule,blockDataIndices,blocks,tapSequence,noiseFD, ...
    noiseVarianceFD,iciNoiseBySymbol,modeName,modulationOrder,predictedBER,targetFeasible,candidateIndex)
N = cfg.system.nFFT; cp = cfg.system.cpLength;
Nt = cfg.system.nTx; Nr = cfg.system.nRx; ns = params.nSymbols;
numData = schedule.numDataSymbols; dataIndices = schedule.dataIndices;
bitsPerQAM = log2(modulationOrder); Xfd = Xpilot;
rng(params.dataSeed+1000+100*candidateIndex,'twister');
if strcmpi(modeName,'spatial-multiplexing')
    txBits = zeros(N*bitsPerQAM,Nt,numData);
    for q = 1:numData
        for tx = 1:Nt
            bits = randi([0 1],N*bitsPerQAM,1);
            txBits(:,tx,q) = bits;
            Xfd(:,tx,dataIndices(q)) = qam_mod(bits,modulationOrder);
        end
    end
else
    numPairs = numData/2;
    txBits = zeros(N*bitsPerQAM,2,numPairs);
    info = complex(zeros(N,2,numPairs));
    for p = 1:numPairs
        for slot = 1:2
            bits = randi([0 1],N*bitsPerQAM,1);
            txBits(:,slot,p) = bits;
            info(:,slot,p) = qam_mod(bits,modulationOrder);
        end
    end
    encoded = alamouti_encode(info);
    for p = 1:numPairs
        Xfd(:,:,dataIndices(2*p-1)) = encoded(:,:,2*p-1);
        Xfd(:,:,dataIndices(2*p)) = encoded(:,:,2*p);
    end
end
rxWaveform = apply_time_varying_mimo_ofdm_channel(ofdm_modulate(Xfd,cp),tapSequence,0);
YfdNoiseless = ofdm_demodulate(rxWaveform,N,cp);
Yfd = YfdNoiseless+noiseFD;
observations = extract_mimo_pilot_observations(Yfd,pilotSymbols,schedule);
[Hest,posteriorVariance] = assemble_burst_channel_estimates_v15_1( ...
    blocks,blockDataIndices,N,Nr,Nt,ns);
rxBits = zeros(size(txBits));
if strcmpi(modeName,'spatial-multiplexing')
    for q = 1:numData
        s = dataIndices(q);
        Hs = reshape(Hest(:,:,:,s),N,Nr,Nt);
        sigmaEff = noiseVarianceFD+sum(posteriorVariance(:,s))+iciNoiseBySymbol(s);
        [Xhat,~] = mimo_equalize(reshape(Yfd(:,:,s),N,Nr,1),Hs,sigmaEff,'mmse');
        for tx = 1:Nt
            rxBits(:,tx,q) = qam_demod(Xhat(:,tx,1),modulationOrder);
        end
    end
else
    for p = 1:numData/2
        s1 = dataIndices(2*p-1); s2 = dataIndices(2*p);
        Ypair = cat(3,Yfd(:,:,s1),Yfd(:,:,s2));
        H1 = reshape(Hest(:,:,:,s1),N,Nr,Nt);
        H2 = reshape(Hest(:,:,:,s2),N,Nr,Nt);
        sigma1 = noiseVarianceFD+sum(posteriorVariance(:,s1))+iciNoiseBySymbol(s1);
        sigma2 = noiseVarianceFD+sum(posteriorVariance(:,s2))+iciNoiseBySymbol(s2);
        [Xhat,~] = alamouti_decode_timevarying(Ypair,H1,H2,[sigma1 sigma2]);
        rxBits(:,1,p) = qam_demod(Xhat(:,1),modulationOrder);
        rxBits(:,2,p) = qam_demod(Xhat(:,2),modulationOrder);
    end
end
bitErrors = nnz(rxBits(:)~=txBits(:)); totalBits = numel(txBits);
[commonMetrics,iciRatio] = measure_candidate_channel_quality_v15_1( ...
    Xfd,YfdNoiseless,Hest,posteriorVariance,tapSequence,N,cp,ns,schedule,observations);
if strcmpi(modeName,'spatial-multiplexing'), streamFactor=Nt; else, streamFactor=1; end
netSE = streamFactor*log2(modulationOrder)*(numData/ns)*N/(N+cp);
result.Mode = char(modeName); result.ModulationOrder=modulationOrder;
result.PredictedBER=predictedBER; result.TargetFeasible=logical(targetFeasible);
result.BER=bitErrors/max(totalBits,1); result.FER=double(bitErrors>0);
result.BitErrors=bitErrors; result.TotalBits=totalBits;
result.NMSE=commonMetrics.NMSE; result.MeanPosteriorVariance=commonMetrics.MeanPosteriorVariance;
result.MeanEffectiveNoiseVariance=mean(noiseVarianceFD+sum(posteriorVariance(:,dataIndices),1)+iciNoiseBySymbol(dataIndices));
result.MeanEstimatedICINoiseVariance=mean(iciNoiseBySymbol(dataIndices));
result.ICI_Ratio=iciRatio; result.PilotContaminationVariance=commonMetrics.PilotContaminationVariance;
result.PilotOverhead=schedule.pilotOverhead; result.NumPilotSymbols=schedule.numPilotSymbols;
result.NumDataSymbols=numData; result.NetSpectralEfficiency_bpsHz=netSE;
result.UncodedGoodput_bpsHz=netSE*(1-result.BER); result.FrameError=bitErrors>0;
if strcmpi(modeName,'spatial-multiplexing')
    result.SpatialMultiplexingFrameRatio=1; result.AlamoutiFrameRatio=0;
else
    result.SpatialMultiplexingFrameRatio=0; result.AlamoutiFrameRatio=1;
end
result.QPSKFrameRatio=double(modulationOrder==4);
result.QAM16FrameRatio=double(modulationOrder==16);
result.TargetFeasibleBlockRatio=double(targetFeasible);
result.CausalReceiverEstimates=true;
end

function [metrics,iciRatio] = measure_candidate_channel_quality_v15_1( ...
    Xfd,YfdNoiseless,Hest,posteriorVariance,tapSequence,N,cp,ns,schedule,observations)
Hdiag=time_varying_ofdm_diagonal_response(tapSequence,N,cp,ns);
dataIndices=schedule.dataIndices;
channelError=Hest(:,:,:,dataIndices)-Hdiag(:,:,:,dataIndices);
selectedTruth=Hdiag(:,:,:,dataIndices);
metrics.NMSE=sum(abs(channelError(:)).^2)/max(sum(abs(selectedTruth(:)).^2),eps);
Ydiag=complex(zeros(size(YfdNoiseless)));
for s=1:ns
    for k=1:N
        Hk=reshape(Hdiag(k,:,: ,s),size(Hdiag,2),size(Hdiag,3));
        xk=reshape(Xfd(k,:,s),size(Xfd,2),1);
        Ydiag(k,:,s)=(Hk*xk).';
    end
end
iciResidual=YfdNoiseless(:,:,dataIndices)-Ydiag(:,:,dataIndices);
iciReference=Ydiag(:,:,dataIndices);
iciRatio=sum(abs(iciResidual(:)).^2)/max(sum(abs(iciReference(:)).^2),eps);
selectedPosterior=posteriorVariance(:,dataIndices);
metrics.MeanPosteriorVariance=mean(selectedPosterior(:));
pilotResidualPower=0; pilotResidualCount=0;
for tx=1:schedule.nTx
    idx=schedule.pilotIndicesByTx{tx};
    for q=1:numel(idx)
        residual=observations{tx}(:,:,q)-Hdiag(:,:,tx,idx(q));
        pilotResidualPower=pilotResidualPower+sum(abs(residual(:)).^2);
        pilotResidualCount=pilotResidualCount+numel(residual);
    end
end
metrics.PilotContaminationVariance=pilotResidualPower/max(pilotResidualCount,1);
end

function [diffuseVariance,meanChannelPower] = channel_variance_config(cfg,modelName)
if strcmpi(char(modelName),'rayleigh')
    diffuseVariance=1; meanChannelPower=1;
elseif strcmpi(char(modelName),'rician')
    K=10^(cfg.channel.ricianK_dB/10);
    diffuseVariance=1/(K+1);
    meanChannelPower=diffuseVariance+cfg.channel.pdp(1)*K/(K+1);
else
    error('simulate_pilot_mapped_adaptive_frame_v15_1:UnsupportedChannel', ...
        'modelName must be rayleigh or rician.');
end
end
