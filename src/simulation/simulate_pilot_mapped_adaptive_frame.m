function frame = simulate_pilot_mapped_adaptive_frame(cfg,modelName,params)
%SIMULATE_PILOT_MAPPED_ADAPTIVE_FRAME Pilot-driven adaptive MIMO-OFDM frame.
% The first orthogonal pilot occasion is used to make a causal frame-level
% decision among {SM, Alamouti} x {QPSK, 16-QAM}. The selected payload is then
% sent through the same sample-varying MIMO FIR waveform channel as V13/V14.
% Later pilots track the channel with causal Jakes-LMMSE for data detection.
%
% The Alamouti receiver uses the estimated channel at both OFDM slots and
% therefore does not require a block-static channel across the codeword.
%
% Required params fields:
%   nSymbols, pilotSpacing, pilotSNRdB, maxDopplerHz, txRho, rxRho,
%   dataSeed, pilotSeed, channelSeed, noiseSeed, targetBER, modulations.
% Optional: estimator (default 'jakes-lmmse').

required = {'nSymbols','pilotSpacing','pilotSNRdB','maxDopplerHz','txRho','rxRho', ...
    'dataSeed','pilotSeed','channelSeed','noiseSeed','targetBER','modulations'};
if ~isstruct(params) || ~all(isfield(params,required))
    error('simulate_pilot_mapped_adaptive_frame:InvalidParams', ...
        'params is missing one or more required fields.');
end
if ~isstruct(cfg) || ~isfield(cfg,'system') || ~isfield(cfg,'channel')
    error('simulate_pilot_mapped_adaptive_frame:InvalidConfig','Invalid project configuration.');
end
if cfg.system.nTx ~= 2 || cfg.system.nRx < 2
    error('simulate_pilot_mapped_adaptive_frame:Topology', ...
        'V15 currently requires two transmit antennas and at least two receive antennas.');
end
if ~isscalar(params.nSymbols) || ~isnumeric(params.nSymbols) || ~isreal(params.nSymbols) || ...
        ~isfinite(params.nSymbols) || params.nSymbols ~= floor(params.nSymbols) || ...
        ~isscalar(params.pilotSpacing) || ~isnumeric(params.pilotSpacing) || ...
        ~isreal(params.pilotSpacing) || ~isfinite(params.pilotSpacing) || ...
        params.pilotSpacing ~= floor(params.pilotSpacing) || ...
        params.nSymbols < 4 || params.pilotSpacing < cfg.system.nTx
    error('simulate_pilot_mapped_adaptive_frame:Schedule', ...
        'nSymbols and pilotSpacing must be valid positive integers.');
end
if ~isscalar(params.pilotSNRdB) || ~isreal(params.pilotSNRdB) || ~isfinite(params.pilotSNRdB) || ...
        ~isscalar(params.maxDopplerHz) || ~isreal(params.maxDopplerHz) || ...
        ~isfinite(params.maxDopplerHz) || params.maxDopplerHz < 0
    error('simulate_pilot_mapped_adaptive_frame:LinkParameters', ...
        'Pilot SNR must be finite and Doppler must be finite and nonnegative.');
end
if ~isscalar(params.txRho) || ~isscalar(params.rxRho) || ...
        ~isreal(params.txRho) || ~isreal(params.rxRho) || ...
        ~isfinite(params.txRho) || ~isfinite(params.rxRho) || ...
        params.txRho < 0 || params.rxRho < 0 || ...
        params.txRho > 0.99 || params.rxRho > 0.99
    error('simulate_pilot_mapped_adaptive_frame:Correlation', ...
        'Spatial-correlation coefficients must be real scalars in [0,0.99].');
end
if ~isscalar(params.targetBER) || ~isreal(params.targetBER) || ...
        ~isfinite(params.targetBER) || params.targetBER <= 0 || params.targetBER >= 0.5
    error('simulate_pilot_mapped_adaptive_frame:TargetBER','targetBER must lie in (0,0.5).');
end
modulations = double(params.modulations(:).');
if isempty(modulations) || ~all(ismember([4 16],modulations)) || any(diff(modulations) <= 0) || ...
        any(modulations < 4) || any(abs(log2(modulations)-round(log2(modulations))) > eps) || ...
        any(abs(log2(sqrt(modulations))-round(log2(sqrt(modulations)))) > eps)
    error('simulate_pilot_mapped_adaptive_frame:Modulations', ...
        'modulations must be ascending square power-of-four QAM orders.');
end
if ~all(isfield(params,{'dataSeed','pilotSeed','channelSeed','noiseSeed'})) || ...
        any(~isfinite([params.dataSeed params.pilotSeed params.channelSeed params.noiseSeed])) || ...
        any([params.dataSeed params.pilotSeed params.channelSeed params.noiseSeed] < 0) || ...
        any([params.dataSeed params.pilotSeed params.channelSeed params.noiseSeed] ~= ...
        floor([params.dataSeed params.pilotSeed params.channelSeed params.noiseSeed]))
    error('simulate_pilot_mapped_adaptive_frame:Seeds','All seeds must be nonnegative integers.');
end
if isfield(params,'estimator')
    estimator = char(params.estimator);
else
    estimator = 'jakes-lmmse';
end
if ~strcmpi(estimator,'jakes-lmmse')
    error('simulate_pilot_mapped_adaptive_frame:Estimator', ...
        'V15 currently uses the causal jakes-lmmse tracker.');
end

N = cfg.system.nFFT;
cp = cfg.system.cpLength;
Nt = cfg.system.nTx;
Nr = cfg.system.nRx;
ns = params.nSymbols;
schedule = mimo_ofdm_pilot_schedule(ns,Nt,params.pilotSpacing);
if mod(schedule.numDataSymbols,2) ~= 0
    error('simulate_pilot_mapped_adaptive_frame:OddDataSlots', ...
        'Alamouti requires an even number of data OFDM symbols after pilot reservation.');
end

cfgFrame = cfg;
cfgFrame.channel.doppler.maxDopplerHz = params.maxDopplerHz;
cfgFrame.channel.spatialCorrelation.enabled = true;
cfgFrame.channel.spatialCorrelation.txRho = params.txRho;
cfgFrame.channel.spatialCorrelation.rxRho = params.rxRho;

% Generate the orthogonal, full-band TDM pilot grid.
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(params.pilotSeed,'twister');
nOccasions = numel(schedule.pilotOccasionStarts);
pilotSymbols = complex(zeros(N,Nt,nOccasions));
Xpilot = complex(zeros(N,Nt,ns));
for tx = 1:Nt
    for q = 1:nOccasions
        knownBits = randi([0 1],2*N,1);
        knownSequence = sqrt(cfg.system.pilotPower)*qam_mod(knownBits,4);
        pilotSymbols(:,tx,q) = knownSequence;
        pilotIndex = schedule.pilotIndicesByTx{tx}(q);
        Xpilot(:,tx,pilotIndex) = knownSequence;
    end
end

blockLength = N+cp;
sampleRateHz = N*cfg.system.subcarrierSpacingHz;
totalSamples = blockLength*ns;
[tapSequence,components,channelMeta] = generate_time_varying_mimo_channel( ...
    cfgFrame,modelName,totalSamples,sampleRateHz,params.channelSeed);

% Match the V13/V14 pilot-SNR convention: noiseVarianceFD is per FFT bin.
if strcmpi(char(modelName),'rayleigh')
    diffuseVariance = 1;
    meanChannelPower = 1;
elseif strcmpi(char(modelName),'rician')
    K = 10^(cfg.channel.ricianK_dB/10);
    diffuseVariance = 1/(K+1);
    meanChannelPower = diffuseVariance + cfg.channel.pdp(1)*K/(K+1);
else
    error('simulate_pilot_mapped_adaptive_frame:UnsupportedChannel', ...
        'modelName must be rayleigh or rician.');
end
noiseVarianceFD = meanChannelPower*10^(-params.pilotSNRdB/10);
rng(params.noiseSeed,'twister');
noiseFD = sqrt(noiseVarianceFD/2)*(randn(N,Nr,ns)+1i*randn(N,Nr,ns));

% The pilot-only waveform gives the first decision occasion. The same pilot
% symbols and noise samples are reused in all full-frame candidate waveforms.
pilotRx = apply_time_varying_mimo_ofdm_channel(ofdm_modulate(Xpilot,cp),tapSequence,0);
pilotY = ofdm_demodulate(pilotRx,N,cp)+noiseFD;
initialObservations = extract_mimo_pilot_observations(pilotY,pilotSymbols,schedule);
Hmean = time_varying_ofdm_diagonal_response(components.los,N,cp,ns);
Tsym = blockLength/(N*cfg.system.subcarrierSpacingHz);
pilotLSNoiseVariance = noiseVarianceFD/cfg.system.pilotPower;
firstDataIndex = schedule.dataIndices(1);
Hdecision = complex(zeros(N,Nr,Nt));
for tx = 1:Nt
    firstPilotIndex = schedule.pilotIndicesByTx{tx}(1);
    obsOne = reshape(initialObservations{tx}(:,:,1),N,Nr,1,1);
    prior4 = reshape(Hmean(:,:,tx,:),N,Nr,1,ns);
    decisionSetup = temporal_channel_tracker_setup(firstPilotIndex,ns, ...
        pilotLSNoiseVariance,params.maxDopplerHz,Tsym,diffuseVariance);
    decisionTrack = temporal_channel_tracking_apply(obsOne,firstPilotIndex,ns, ...
        'jakes-lmmse',prior4,decisionSetup);
    Hdecision(:,:,tx) = decisionTrack(:,:,1,firstDataIndex);
end

% Score each MIMO/QAM candidate using estimated post-detection SINR and the
% actual pilot and CP overhead of this frame's schedule.
streamFactors = [Nt;1];
spectralEfficiency = zeros(2,numel(modulations));
for modeIndex = 1:2
    for m = 1:numel(modulations)
        spectralEfficiency(modeIndex,m) = streamFactors(modeIndex)*log2(modulations(m))* ...
            (schedule.numDataSymbols/ns)*N/(N+cp);
    end
end
prediction = predict_pilot_mapped_mcs(Hdecision,noiseVarianceFD, ...
    modulations,params.targetBER,spectralEfficiency);
choice = prediction.choice;

% Run the selected configuration and two fixed reference links using the same
% channel, pilots and thermal-noise realization (paired comparison).
adaptiveResult = simulate_one_candidate(cfgFrame,modelName,params, ...
    Xpilot,pilotSymbols,schedule,tapSequence,noiseFD,Hmean,noiseVarianceFD, ...
    choice.Mode,choice.ModulationOrder,choice.PredictedBER,choice.TargetFeasible,1);
fixedSM16 = simulate_one_candidate(cfgFrame,modelName,params, ...
    Xpilot,pilotSymbols,schedule,tapSequence,noiseFD,Hmean,noiseVarianceFD, ...
    'spatial-multiplexing',16,prediction.predictedBER(1,find(modulations==16,1)), ...
    prediction.predictedBER(1,find(modulations==16,1))<=params.targetBER,2);
fixedAlamoutiQPSK = simulate_one_candidate(cfgFrame,modelName,params, ...
    Xpilot,pilotSymbols,schedule,tapSequence,noiseFD,Hmean,noiseVarianceFD, ...
    'alamouti-stbc',4,prediction.predictedBER(2,find(modulations==4,1)), ...
    prediction.predictedBER(2,find(modulations==4,1))<=params.targetBER,3);

frame.schedule = schedule;
frame.channelMeta = channelMeta;
frame.noiseVarianceFD = noiseVarianceFD;
frame.targetBER = params.targetBER;
frame.estimator = estimator;
frame.Hdecision = Hdecision;
frame.prediction = prediction;
frame.choice = choice;
frame.adaptive = adaptiveResult;
frame.fixedSM16QAM = fixedSM16;
frame.fixedAlamoutiQPSK = fixedAlamoutiQPSK;
end

function result = simulate_one_candidate(cfg,modelName,params,Xpilot,pilotSymbols, ...
    schedule,tapSequence,noiseFD,Hmean,noiseVarianceFD,modeName,modulationOrder, ...
    predictedBER,targetFeasible,candidateIndex)
N = cfg.system.nFFT;
cp = cfg.system.cpLength;
Nt = cfg.system.nTx;
Nr = cfg.system.nRx;
ns = params.nSymbols;
numData = schedule.numDataSymbols;
dataIndices = schedule.dataIndices;
bitsPerQAM = log2(modulationOrder);
Xfd = Xpilot;
rng(params.dataSeed+100*candidateIndex,'twister');

if strcmpi(modeName,'spatial-multiplexing')
    streamCount = Nt;
    txBits = zeros(N*bitsPerQAM,Nt,numData);
    for q = 1:numData
        symbolIndex = dataIndices(q);
        for tx = 1:Nt
            bits = randi([0 1],N*bitsPerQAM,1);
            txBits(:,tx,q) = bits;
            Xfd(:,tx,symbolIndex) = qam_mod(bits,modulationOrder);
        end
    end
elseif strcmpi(modeName,'alamouti-stbc')
    streamCount = 1;
    numPairs = numData/2;
    txBits = zeros(N*bitsPerQAM,2,numPairs);
    informationGrid = complex(zeros(N,2,numPairs));
    for pairIndex = 1:numPairs
        for informationSlot = 1:2
            bits = randi([0 1],N*bitsPerQAM,1);
            txBits(:,informationSlot,pairIndex) = bits;
            informationGrid(:,informationSlot,pairIndex) = qam_mod(bits,modulationOrder);
        end
    end
    encodedGrid = alamouti_encode(informationGrid);
    for pairIndex = 1:numPairs
        firstDataPosition = 2*pairIndex-1;
        secondDataPosition = 2*pairIndex;
        Xfd(:,:,dataIndices(firstDataPosition)) = encodedGrid(:,:,2*pairIndex-1);
        Xfd(:,:,dataIndices(secondDataPosition)) = encodedGrid(:,:,2*pairIndex);
    end
else
    error('simulate_pilot_mapped_adaptive_frame:UnknownMode','Unknown candidate mode.');
end

rxWaveform = apply_time_varying_mimo_ofdm_channel(ofdm_modulate(Xfd,cp),tapSequence,0);
YfdNoiseless = ofdm_demodulate(rxWaveform,N,cp);
Yfd = YfdNoiseless+noiseFD;
observationsByTx = extract_mimo_pilot_observations(Yfd,pilotSymbols,schedule);
Hest = complex(zeros(N,Nr,Nt,ns));
posteriorVariance = nan(Nt,ns);
Tsym = (N+cp)/(N*cfg.system.subcarrierSpacingHz);
diffuseVariance = get_diffuse_variance(cfg,modelName);
for tx = 1:Nt
    pilotIndices = schedule.pilotIndicesByTx{tx};
    obs4 = reshape(observationsByTx{tx},N,Nr,1,numel(pilotIndices));
    prior4 = reshape(Hmean(:,:,tx,:),N,Nr,1,ns);
    setup = temporal_channel_tracker_setup(pilotIndices,ns, ...
        noiseVarianceFD/cfg.system.pilotPower,params.maxDopplerHz,Tsym,diffuseVariance);
    [localEstimate,localPosterior] = temporal_channel_tracking_apply( ...
        obs4,pilotIndices,ns,'jakes-lmmse',prior4,setup);
    Hest(:,:,tx,:) = localEstimate;
    posteriorVariance(tx,:) = localPosterior(:).';
end

% Per-mode hard decoding.
rxBits = zeros(size(txBits));
if strcmpi(modeName,'spatial-multiplexing')
    for q = 1:numData
        symbolIndex = dataIndices(q);
        Hsymbol = reshape(Hest(:,:,:,symbolIndex),N,Nr,Nt);
        [Xhat,~] = mimo_equalize(reshape(Yfd(:,:,symbolIndex),N,Nr,1), ...
            Hsymbol,noiseVarianceFD,'mmse');
        for tx = 1:Nt
            rxBits(:,tx,q) = qam_demod(Xhat(:,tx,1),modulationOrder);
        end
    end
else
    numPairs = numData/2;
    for pairIndex = 1:numPairs
        firstDataPosition = 2*pairIndex-1;
        secondDataPosition = 2*pairIndex;
        firstSymbol = dataIndices(firstDataPosition);
        secondSymbol = dataIndices(secondDataPosition);
        Ypair = complex(zeros(N,Nr,2));
        Ypair(:,:,1) = Yfd(:,:,firstSymbol);
        Ypair(:,:,2) = Yfd(:,:,secondSymbol);
        Hfirst = reshape(Hest(:,:,:,firstSymbol),N,Nr,Nt);
        Hsecond = reshape(Hest(:,:,:,secondSymbol),N,Nr,Nt);
        [Xhat,~] = alamouti_decode_timevarying(Ypair,Hfirst,Hsecond,noiseVarianceFD);
        rxBits(:,1,pairIndex) = qam_demod(Xhat(:,1),modulationOrder);
        rxBits(:,2,pairIndex) = qam_demod(Xhat(:,2),modulationOrder);
    end
end

bitErrors = nnz(rxBits(:)~=txBits(:));
totalBits = numel(txBits);
BER = bitErrors/totalBits;

Hdiag = time_varying_ofdm_diagonal_response(tapSequence,N,cp,ns);
selectedTruth = Hdiag(:,:,:,dataIndices);
selectedEstimate = Hest(:,:,:,dataIndices);
channelError = selectedEstimate-selectedTruth;
NMSE = sum(abs(channelError(:)).^2)/max(sum(abs(selectedTruth(:)).^2),eps);

Ydiag = complex(zeros(size(YfdNoiseless)));
for s = 1:ns
    for k = 1:N
        Hk = reshape(Hdiag(k,:,:,s),Nr,Nt);
        xk = reshape(Xfd(k,:,s),Nt,1);
        Ydiag(k,:,s) = (Hk*xk).';
    end
end
iciResidual = YfdNoiseless(:,:,dataIndices)-Ydiag(:,:,dataIndices);
iciReference = Ydiag(:,:,dataIndices);
iciRatio = sum(abs(iciResidual(:)).^2)/max(sum(abs(iciReference(:)).^2),eps);

pilotContaminationSum = 0;
pilotContaminationCount = 0;
for tx = 1:Nt
    idx = schedule.pilotIndicesByTx{tx};
    for q = 1:numel(idx)
        residual = observationsByTx{tx}(:,:,q)-Hdiag(:,:,tx,idx(q));
        pilotContaminationSum = pilotContaminationSum+sum(abs(residual(:)).^2);
        pilotContaminationCount = pilotContaminationCount+numel(residual);
    end
end
pilotContaminationVariance = pilotContaminationSum/max(pilotContaminationCount,1);

if streamCount == Nt
    streamFactor = Nt;
else
    streamFactor = 1;
end
netSE = streamFactor*log2(modulationOrder)*(numData/ns)*N/(N+cp);
result.Mode = char(modeName);
result.ModulationOrder = modulationOrder;
result.PredictedBER = predictedBER;
result.TargetFeasible = logical(targetFeasible);
result.BER = BER;
result.FER = double(bitErrors>0);
result.BitErrors = bitErrors;
result.TotalBits = totalBits;
result.NMSE = NMSE;
posteriorSelection = posteriorVariance(:,dataIndices(:));
result.MeanPosteriorVariance = mean(posteriorSelection(:));
result.ICI_Ratio = iciRatio;
result.PilotContaminationVariance = pilotContaminationVariance;
result.PilotOverhead = schedule.pilotOverhead;
result.NumPilotSymbols = schedule.numPilotSymbols;
result.NumDataSymbols = numData;
result.NetSpectralEfficiency_bpsHz = netSE;
result.UncodedGoodput_bpsHz = netSE*(1-BER);
result.FrameError = bitErrors>0;
end

function diffuseVariance = get_diffuse_variance(cfg,modelName)
if strcmpi(char(modelName),'rayleigh')
    diffuseVariance = 1;
elseif strcmpi(char(modelName),'rician')
    K = 10^(cfg.channel.ricianK_dB/10);
    diffuseVariance = 1/(K+1);
else
    error('simulate_pilot_mapped_adaptive_frame:UnsupportedChannel', ...
        'modelName must be rayleigh or rician.');
end
end
