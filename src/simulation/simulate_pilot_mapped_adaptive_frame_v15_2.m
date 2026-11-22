function frame = simulate_pilot_mapped_adaptive_frame_v15_2(cfg,modelName,params)
%SIMULATE_PILOT_MAPPED_ADAPTIVE_FRAME_V15_2 Block-adaptive packet-coded MIMO-OFDM.
% Inherits V15.1 improvements:
%  1) A fresh MIMO-mode/QAM decision is made after each orthogonal pilot occasion.
%  2) The decision uses only pilots observed by that occasion (strictly causal).
%  3) Per-symbol temporal-LMMSE posterior variance is added to effective noise.
%  4) Residual ICI is estimated from causal pilot innovations, with a small-
%     Doppler analytical floor at the first pilot occasion.
%  5) SM and generalized two-slot Alamouti predictions use the predicted
%     channel at the actual data symbols in the upcoming burst.
%
% Both transmission strategies are evaluated with uncoded+CRC and rate-1/2
% convolutional-FEC+CRC packets. All candidates share the channel, pilot grid,
% causal CSI estimates, and receiver thermal-noise realization.

required = {'nSymbols','pilotSpacing','pilotSNRdB','maxDopplerHz','txRho','rxRho', ...
    'dataSeed','pilotSeed','channelSeed','noiseSeed','targetBER','modulations'};
if ~isstruct(params) || ~all(isfield(params,required))
    error('simulate_pilot_mapped_adaptive_frame_v15_2:InvalidParams', ...
        'params is missing one or more required fields.');
end
if ~isstruct(cfg) || ~isfield(cfg,'system') || ~isfield(cfg,'channel')
    error('simulate_pilot_mapped_adaptive_frame_v15_2:InvalidConfig','Invalid project configuration.');
end
if cfg.system.nTx ~= 2 || cfg.system.nRx < 2
    error('simulate_pilot_mapped_adaptive_frame_v15_2:Topology', ...
        'V15.2 requires two transmit antennas and at least two receive antennas.');
end
if ~isscalar(params.nSymbols) || ~isnumeric(params.nSymbols) || ~isreal(params.nSymbols) || ...
        ~isfinite(params.nSymbols) || params.nSymbols ~= floor(params.nSymbols) || ...
        ~isscalar(params.pilotSpacing) || ~isnumeric(params.pilotSpacing) || ...
        ~isreal(params.pilotSpacing) || ~isfinite(params.pilotSpacing) || ...
        params.pilotSpacing ~= floor(params.pilotSpacing) || ...
        params.nSymbols < 4 || params.pilotSpacing < cfg.system.nTx
    error('simulate_pilot_mapped_adaptive_frame_v15_2:Schedule', ...
        'nSymbols and pilotSpacing must be valid integers.');
end
if ~isscalar(params.pilotSNRdB) || ~isreal(params.pilotSNRdB) || ~isfinite(params.pilotSNRdB) || ...
        ~isscalar(params.maxDopplerHz) || ~isreal(params.maxDopplerHz) || ...
        ~isfinite(params.maxDopplerHz) || params.maxDopplerHz < 0
    error('simulate_pilot_mapped_adaptive_frame_v15_2:LinkParameters', ...
        'Pilot SNR must be finite and Doppler must be nonnegative.');
end
if ~isscalar(params.txRho) || ~isscalar(params.rxRho) || ...
        ~isreal(params.txRho) || ~isreal(params.rxRho) || ...
        ~isfinite(params.txRho) || ~isfinite(params.rxRho) || ...
        params.txRho < 0 || params.rxRho < 0 || params.txRho > 0.99 || params.rxRho > 0.99
    error('simulate_pilot_mapped_adaptive_frame_v15_2:Correlation', ...
        'Spatial-correlation coefficients must lie in [0,0.99].');
end
if ~isscalar(params.targetBER) || ~isreal(params.targetBER) || ...
        ~isfinite(params.targetBER) || params.targetBER <= 0 || params.targetBER >= 0.5
    error('simulate_pilot_mapped_adaptive_frame_v15_2:TargetBER', ...
        'targetBER must lie in (0,0.5).');
end
modulations = double(params.modulations(:).');
if isempty(modulations) || ~all(ismember([4 16],modulations)) || ...
        any(diff(modulations) <= 0) || ...
        any(abs(log2(modulations)-round(log2(modulations))) > eps) || ...
        any(abs(log2(sqrt(modulations))-round(log2(sqrt(modulations)))) > eps)
    error('simulate_pilot_mapped_adaptive_frame_v15_2:Modulations', ...
        'modulations must be ascending square power-of-four QAM orders including 4 and 16.');
end
seeds = [params.dataSeed params.pilotSeed params.channelSeed params.noiseSeed];
if any(~isfinite(seeds)) || any(seeds < 0) || any(seeds ~= floor(seeds))
    error('simulate_pilot_mapped_adaptive_frame_v15_2:Seeds','All seeds must be nonnegative integers.');
end
if isfield(params,'estimator') && ~strcmpi(char(params.estimator),'jakes-lmmse')
    error('simulate_pilot_mapped_adaptive_frame_v15_2:Estimator', ...
        'V15.2 currently uses the causal Jakes-LMMSE tracker.');
end

N = cfg.system.nFFT; cp = cfg.system.cpLength;
Nt = cfg.system.nTx; Nr = cfg.system.nRx; ns = params.nSymbols;
schedule = mimo_ofdm_pilot_schedule(ns,Nt,params.pilotSpacing);
blockDataIndices = build_data_blocks_v15_2(schedule,Nt,ns);

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
    [Hprovisional,~] = track_available_pilots_v15_2( ...
        observationsByTx,schedule,b,Hmean,provisionalHistory, ...
        cfg,params,ns,diffuseVariance,Tsym,dataIdx,noiseVarianceFD);
    [iciVariance,iciDetail] = estimate_pilot_ici_variance( ...
        observationsByTx,schedule,b,Hmean,noiseVarianceFD, ...
        cfg.system.pilotPower,params.maxDopplerHz,Tsym,diffuseVariance, ...
        cfg.system.subcarrierSpacingHz,Hprovisional,previousHistory);
    iciNoiseByOccasion(b) = iciVariance;
    [Hburst,posteriorByTx] = track_available_pilots_v15_2( ...
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

% Simulate each physical strategy in uncoded-CRC and convolutional-CRC modes.
strategyNames = {'adaptive','fixed-sm-16qam','fixed-alamouti-qpsk'};
codingModes = {'uncoded-crc','conv-crc-1/2'};
sm16Index = find(modulations==16,1);
qpskIndex = find(modulations==4,1);
sm16Pred = mean(arrayfun(@(x)x.prediction.predictedBER(1,sm16Index),blocks));
alQpskPred = mean(arrayfun(@(x)x.prediction.predictedBER(2,qpskIndex),blocks));
frame.candidates = cell(numel(strategyNames)*numel(codingModes),1);
for si = 1:numel(strategyNames)
    switch si
        case 1
            modeByBlock = arrayfun(@(x)char(x.choice.Mode),blocks,'UniformOutput',false);
            modByBlock = arrayfun(@(x)x.choice.ModulationOrder,blocks);
            predByBlock = arrayfun(@(x)x.choice.PredictedBER,blocks);
            feasByBlock = arrayfun(@(x)x.choice.TargetFeasible,blocks);
        case 2
            modeByBlock = repmat({'spatial-multiplexing'},numel(blocks),1);
            modByBlock = 16*ones(1,numel(blocks));
            predByBlock = sm16Pred*ones(1,numel(blocks));
            feasByBlock = predByBlock<=params.targetBER;
        case 3
            modeByBlock = repmat({'alamouti-stbc'},numel(blocks),1);
            modByBlock = 4*ones(1,numel(blocks));
            predByBlock = alQpskPred*ones(1,numel(blocks));
            feasByBlock = predByBlock<=params.targetBER;
    end
    modeByBlock=modeByBlock(:).';
    modByBlock=modByBlock(:).';
    predByBlock=predByBlock(:).';
    feasByBlock=logical(feasByBlock(:).');
    for ci = 1:numel(codingModes)
        candidate = simulate_candidate_v15_2(cfgFrame,params, ...
            Xpilot,pilotSymbols,schedule,blockDataIndices,blocks,tapSequence,noiseFD, ...
            noiseVarianceFD,iciNoiseBySymbol,modeByBlock,modByBlock,predByBlock, ...
            feasByBlock,codingModes{ci},1);
        frame.candidates{(si-1)*numel(codingModes)+ci} = candidate;
    end
end

frame.schedule = schedule;
frame.dataBlocks = blockDataIndices;
frame.channelMeta = channelMeta;
frame.noiseVarianceFD = noiseVarianceFD;
frame.targetBER = params.targetBER;
frame.estimator = 'jakes-lmmse';
frame.blocks = blocks;
frame.choice = blocks(1).choice;
frame.prediction = blocks(1).prediction;
frame.strategyNames = strategyNames;
frame.codingModes = codingModes;
end

function [Hburst,posteriorByTx] = track_available_pilots_v15_2( ...
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

function blocks = build_data_blocks_v15_2(schedule,nTx,nSymbols)
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
        error('simulate_pilot_mapped_adaptive_frame_v15_2:OddDataSlots', ...
            'Each pilot-to-pilot data burst must contain an even number of OFDM symbols.');
    end
    blocks{b} = idx;
end
end

function result = simulate_candidate_v15_2(cfg,params,Xpilot,pilotSymbols,schedule, ...
    blockDataIndices,blocks,tapSequence,noiseFD,noiseVarianceFD,iciNoiseBySymbol, ...
    modeByBlock,modByBlock,predByBlock,feasByBlock,codingMode,candidateSeed)
N=cfg.system.nFFT; cp=cfg.system.cpLength; Nt=cfg.system.nTx;
Nr=cfg.system.nRx; ns=params.nSymbols; Xfd=Xpilot;
numBlocks=numel(blockDataIndices);
packetMeta=cell(numBlocks,1); payloadBitsByBlock=cell(numBlocks,1);
grossBitsByBlock=zeros(numBlocks,1); blockGrossSE=zeros(numBlocks,1);
rng(params.dataSeed+candidateSeed,'twister');

% Packetize every pilot-to-pilot data burst. One packet spans the active
% spatial layers of that burst and carries a real CRC; FEC mode uses a
% terminated (7,5) rate-1/2 convolutional code.
for b=1:numBlocks
    idx=blockDataIndices{b}; M=modByBlock(b); bps=log2(M);
    if strcmpi(modeByBlock{b},'spatial-multiplexing')
        capacity=N*bps*Nt*numel(idx);
        packet=packet_crc_encode_v15_2(capacity,codingMode, ...
            params.dataSeed+candidateSeed*1000+b);
        bitGrid=reshape(packet.EncodedBits,N*bps,Nt,numel(idx));
        for q=1:numel(idx)
            for tx=1:Nt
                Xfd(:,tx,idx(q))=qam_mod(bitGrid(:,tx,q),M);
            end
        end
        blockGrossSE(b)=Nt*bps;
    else
        numPairs=numel(idx)/2;
        if mod(numel(idx),2)~=0
            error('simulate_pilot_mapped_adaptive_frame_v15_2:AlamoutiPairs', ...
                'Alamouti data bursts must have an even number of OFDM symbols.');
        end
        capacity=N*bps*numel(idx);
        packet=packet_crc_encode_v15_2(capacity,codingMode, ...
            params.dataSeed+candidateSeed*1000+b);
        bitGrid=reshape(packet.EncodedBits,N*bps,2,numPairs);
        info=complex(zeros(N,2,numPairs));
        for pair=1:numPairs
            for slot=1:2
                info(:,slot,pair)=qam_mod(bitGrid(:,slot,pair),M);
            end
        end
        encoded=alamouti_encode(info);
        for pair=1:numPairs
            Xfd(:,:,idx(2*pair-1))=encoded(:,:,2*pair-1);
            Xfd(:,:,idx(2*pair))=encoded(:,:,2*pair);
        end
        blockGrossSE(b)=bps;
    end
    packetMeta{b}=packet;
    payloadBitsByBlock{b}=packet.PayloadBits;
    grossBitsByBlock(b)=packet.EncodedCapacityBits;
end

rxWaveform=apply_time_varying_mimo_ofdm_channel(ofdm_modulate(Xfd,cp),tapSequence,0);
YfdNoiseless=ofdm_demodulate(rxWaveform,N,cp);
Yfd=YfdNoiseless+noiseFD;
observations=extract_mimo_pilot_observations(Yfd,pilotSymbols,schedule);
[Hest,posteriorVariance]=assemble_burst_channel_estimates_v15_1( ...
    blocks,blockDataIndices,N,Nr,Nt,ns);

payloadErrors=0; totalPayloadBits=0; codedBitErrors=0; totalCodedBits=0;
packetErrors=0; undetectedErrors=0; acceptedPayloadBits=0;
for b=1:numBlocks
    idx=blockDataIndices{b}; modeName=modeByBlock{b}; M=modByBlock(b);
    bps=log2(M);
    if strcmpi(modeName,'spatial-multiplexing')
        rxLlrGrid=zeros(N*bps,Nt,numel(idx));
        for q=1:numel(idx)
            s=idx(q); Hs=reshape(Hest(:,:,:,s),N,Nr,Nt);
            sigmaEff=noiseVarianceFD+sum(posteriorVariance(:,s))+iciNoiseBySymbol(s);
            [Xhat,postSinr]=mimo_equalize(reshape(Yfd(:,:,s),N,Nr,1),Hs,sigmaEff,'mmse');
            for tx=1:Nt
                symbolVariance=1./max(postSinr(:,tx,1),1e-12);
                rxLlrGrid(:,tx,q)=qam_soft_demod_llr(Xhat(:,tx,1),M,symbolVariance);
            end
        end
    else
        numPairs=numel(idx)/2;
        rxLlrGrid=zeros(N*bps,2,numPairs);
        for pair=1:numPairs
            s1=idx(2*pair-1); s2=idx(2*pair);
            Ypair=cat(3,Yfd(:,:,s1),Yfd(:,:,s2));
            H1=reshape(Hest(:,:,:,s1),N,Nr,Nt);
            H2=reshape(Hest(:,:,:,s2),N,Nr,Nt);
            sigma1=noiseVarianceFD+sum(posteriorVariance(:,s1))+iciNoiseBySymbol(s1);
            sigma2=noiseVarianceFD+sum(posteriorVariance(:,s2))+iciNoiseBySymbol(s2);
            [Xhat,postSinr]=alamouti_decode_timevarying(Ypair,H1,H2,[sigma1 sigma2]);
            for layer=1:2
                symbolVariance=1./max(postSinr(:,layer),1e-12);
                rxLlrGrid(:,layer,pair)=qam_soft_demod_llr(Xhat(:,layer),M,symbolVariance);
            end
        end
    end
    decoded=packet_soft_decode_v15_2(rxLlrGrid(:),packetMeta{b});
    payloadErrors=payloadErrors+decoded.PayloadErrors;
    totalPayloadBits=totalPayloadBits+decoded.PayloadLength;
    codedBitErrors=codedBitErrors+decoded.CodedBitErrors;
    totalCodedBits=totalCodedBits+decoded.ReceivedEncodedBits;
    packetErrors=packetErrors+double(decoded.PacketError);
    undetectedErrors=undetectedErrors+double(decoded.UndetectedPayloadErrors);
    if decoded.CRCValid
        acceptedPayloadBits=acceptedPayloadBits+decoded.PayloadLength;
    end
end

[commonMetrics,iciRatio]=measure_candidate_channel_quality_v15_2( ...
    Xfd,YfdNoiseless,Hest,posteriorVariance,tapSequence,N,cp,ns,schedule,observations);
rawGrossSE=sum(blockGrossSE.*cellfun(@numel,blockDataIndices(:)))/ns*N/(N+cp);
payloadSE=totalPayloadBits/(N*ns)*N/(N+cp);
crcGoodput=acceptedPayloadBits/(N*ns)*N/(N+cp);
result.Mode='block-adaptive';
if numel(unique(modeByBlock))==1, result.Mode=modeByBlock{1}; end
result.ModulationOrder=NaN; if numel(unique(modByBlock))==1, result.ModulationOrder=modByBlock(1); end
result.ModeByBlock=modeByBlock; result.ModulationOrderByBlock=modByBlock;
slotWeights=cellfun(@numel,blockDataIndices(:)); slotWeights=slotWeights(:).';
result.SpatialMultiplexingFrameRatio=sum(strcmp(modeByBlock,'spatial-multiplexing').*slotWeights)/schedule.numDataSymbols;
result.AlamoutiFrameRatio=1-result.SpatialMultiplexingFrameRatio;
result.QPSKFrameRatio=sum((modByBlock==4).*slotWeights)/schedule.numDataSymbols;
result.QAM16FrameRatio=sum((modByBlock==16).*slotWeights)/schedule.numDataSymbols;
result.PredictedBER=sum(predByBlock(:).'.*slotWeights)/sum(slotWeights);
result.TargetFeasible=all(feasByBlock);
result.TargetFeasibleBlockRatio=sum(feasByBlock(:).'.*slotWeights)/schedule.numDataSymbols;
result.BER=payloadErrors/max(totalPayloadBits,1);
result.PayloadBER=result.BER;
result.CodedBER=codedBitErrors/max(totalCodedBits,1);
result.BitErrors=payloadErrors; result.TotalBits=totalPayloadBits;
result.CodedBitErrors=codedBitErrors; result.TotalCodedBits=totalCodedBits;
result.PacketErrorRate=packetErrors/max(numBlocks,1);
result.PacketErrors=packetErrors; result.PacketCount=numBlocks;
result.UndetectedPacketErrors=undetectedErrors;
result.FER=double(packetErrors>0); result.FrameError=packetErrors>0;
result.CRCAcceptedPayloadBits=acceptedPayloadBits;
result.PayloadBits=totalPayloadBits;
result.NMSE=commonMetrics.NMSE;
result.MeanPosteriorVariance=commonMetrics.MeanPosteriorVariance;
result.MeanEffectiveNoiseVariance=mean(noiseVarianceFD+sum(posteriorVariance(:,schedule.dataIndices),1)+iciNoiseBySymbol(schedule.dataIndices));
result.MeanEstimatedICINoiseVariance=mean(iciNoiseBySymbol(schedule.dataIndices));
result.ICI_Ratio=iciRatio; result.PilotContaminationVariance=commonMetrics.PilotContaminationVariance;
result.PilotOverhead=schedule.pilotOverhead; result.NumPilotSymbols=schedule.numPilotSymbols;
result.NumDataSymbols=schedule.numDataSymbols;
result.GrossSpectralEfficiency_bpsHz=rawGrossSE;
result.NetSpectralEfficiency_bpsHz=payloadSE;
result.BitGoodput_bpsHz=payloadSE*(1-result.BER); % expected correct payload-bit rate
result.UncodedGoodput_bpsHz=result.BitGoodput_bpsHz; % legacy compatibility alias
result.CRCGoodput_bpsHz=crcGoodput; % only CRC-accepted complete packets count
result.NumAdaptationBlocks=numBlocks; result.CausalReceiverEstimates=true;
result.CodingMode=char(codingMode);
result.QAMDemapper='max-log soft LLR';
if strcmpi(codingMode,'conv-crc-1/2')
    result.FECDecoder='max-log soft-decision Viterbi';
else
    result.FECDecoder='not used (uncoded packet)';
end
result.CodeRate=packetMeta{1}.CodeRate;
end

function [metrics,iciRatio] = measure_candidate_channel_quality_v15_2( ...
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
    error('simulate_pilot_mapped_adaptive_frame_v15_2:UnsupportedChannel', ...
        'modelName must be rayleigh or rician.');
end
end
