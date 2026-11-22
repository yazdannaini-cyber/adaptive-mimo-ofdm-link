function frame = simulate_pilot_mapped_ofdm_frame(cfg,modelName,params)
%SIMULATE_PILOT_MAPPED_OFDM_FRAME End-to-end OFDM frame with mapped MIMO pilots.
% This experiment uses actual time-domain OFDM modulation, a sample-varying
% MIMO FIR channel, CP removal/FFT, pilot extraction, temporal tracking, and
% per-subcarrier detection for estimated/diagonal-CSI methods, plus a full-band
% ICI-aware oracle detector using the exact time-varying tap sequence.
%
% Required params fields:
%   nSymbols, pilotSpacing, modulationOrder, pilotSNRdB, modelName-independent
%   maxDopplerHz, txRho, rxRho, dataSeed, pilotSeed, channelSeed, noiseSeed.

required = {'nSymbols','pilotSpacing','modulationOrder','pilotSNRdB', ...
    'maxDopplerHz','txRho','rxRho','dataSeed','pilotSeed','channelSeed','noiseSeed'};
if ~isstruct(params) || ~all(isfield(params,required))
    error('simulate_pilot_mapped_ofdm_frame:InvalidParams', ...
        'params is missing one or more required fields.');
end
if ~isstruct(cfg) || ~isfield(cfg,'system') || ~isfield(cfg,'channel')
    error('simulate_pilot_mapped_ofdm_frame:InvalidConfig','Invalid project configuration.');
end

N = cfg.system.nFFT;
cp = cfg.system.cpLength;
Nt = cfg.system.nTx;
Nr = cfg.system.nRx;
ns = params.nSymbols;
if Nt ~= 2
    error('simulate_pilot_mapped_ofdm_frame:Topology', ...
        'Pilot-mapped time-division pilot mapping currently targets the project 2x2 MIMO topology.');
end

cfgFrame = cfg;
cfgFrame.channel.doppler.maxDopplerHz = params.maxDopplerHz;
cfgFrame.channel.spatialCorrelation.enabled = true;
cfgFrame.channel.spatialCorrelation.txRho = params.txRho;
cfgFrame.channel.spatialCorrelation.rxRho = params.rxRho;

[Xfd,dataBits,pilotSymbols,schedule] = build_mimo_ofdm_pilot_data_grid( ...
    N,Nt,ns,params.pilotSpacing,params.modulationOrder, ...
    cfg.system.pilotPower,params.dataSeed,params.pilotSeed);
txWaveform = ofdm_modulate(Xfd,cp);
blockLength = N+cp;
sampleRateHz = N*cfg.system.subcarrierSpacingHz;
totalSamples = blockLength*ns;
[tapSequence,components,channelMeta] = generate_time_varying_mimo_channel( ...
    cfgFrame,modelName,totalSamples,sampleRateHz,params.channelSeed);
rxWaveformNoiseless = apply_time_varying_mimo_ofdm_channel(txWaveform,tapSequence,0);
YfdNoiseless = ofdm_demodulate(rxWaveformNoiseless,N,cp);
Hdiag = time_varying_ofdm_diagonal_response(tapSequence,N,cp,ns);

% The FFT-domain noise variance equals N times the time-sample variance.
% Adding white Gaussian noise after the FFT is statistically equivalent for
% a CP-removed unitary-free N-point FFT, and keeps the requested per-tone SNR explicit.
if strcmpi(char(modelName),'rayleigh')
    diffuseVariance = 1;
    meanChannelPower = 1;
elseif strcmpi(char(modelName),'rician')
    K = 10^(cfg.channel.ricianK_dB/10);
    diffuseVariance = 1/(K+1);
    meanChannelPower = diffuseVariance + cfg.channel.pdp(1)*K/(K+1);
else
    error('simulate_pilot_mapped_ofdm_frame:UnsupportedChannel', ...
        'modelName must be rayleigh or rician.');
end
noiseVarianceFD = meanChannelPower*10^(-params.pilotSNRdB/10);
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(params.noiseSeed,'twister');
noiseFD = sqrt(noiseVarianceFD/2)*( ...
    randn(size(YfdNoiseless))+1i*randn(size(YfdNoiseless)));
Yfd = YfdNoiseless+noiseFD;

observationsByTx = extract_mimo_pilot_observations(Yfd,pilotSymbols,schedule);
Hmean = time_varying_ofdm_diagonal_response(components.los,N,cp,ns);
Tsym = blockLength/(N*cfg.system.subcarrierSpacingHz);

% Build the diagonal frequency-domain prediction used to quantify ICI.
Ydiag = complex(zeros(size(YfdNoiseless)));
for s = 1:ns
    for k = 1:N
        Hk = reshape(Hdiag(k,:,:,s),Nr,Nt);
        xk = reshape(Xfd(k,:,s),Nt,1);
        Ydiag(k,:,s) = (Hk*xk).';
    end
end
iciResidual = YfdNoiseless(:,:,schedule.dataIndices)-Ydiag(:,:,schedule.dataIndices);
iciReference = Ydiag(:,:,schedule.dataIndices);
iciRatio = sum(abs(iciResidual(:)).^2)/max(sum(abs(iciReference(:)).^2),eps);

% Pilot LS contamination includes thermal noise plus intra-symbol ICI.
pilotContaminationSum = 0;
pilotContaminationCount = 0;
for tx = 1:Nt
    idx = schedule.pilotIndicesByTx{tx};
    obs = observationsByTx{tx};
    for q = 1:numel(idx)
        residual = obs(:,:,q)-Hdiag(:,:,tx,idx(q));
        pilotContaminationSum = pilotContaminationSum+sum(abs(residual(:)).^2);
        pilotContaminationCount = pilotContaminationCount+numel(residual);
    end
end
pilotContaminationVariance = pilotContaminationSum/max(pilotContaminationCount,1);

methodNames = {'hold-last','linear','jakes-lmmse','perfect-csi','perfect-csi-ici'};
methods = repmat(struct('name','','Hest',[],'NMSE',NaN,'MSEPerCoefficient',NaN, ...
    'PosteriorMSEPerCoefficient',NaN,'PosteriorCalibrationRatio',NaN, ...
    'BitErrors',0,'TotalBits',numel(dataBits),'BER',NaN,'FrameError',false), ...
    numel(methodNames),1);
posteriorMaps = nan(N,Nr,Nt,ns);

for im = 1:numel(methodNames)
    method = methodNames{im};
    if strcmp(method,'perfect-csi') || strcmp(method,'perfect-csi-ici')
        % Both oracle references know the exact diagonal channel. The ICI-aware
        % method additionally uses the full off-diagonal frequency-coupling matrix.
        Hest = Hdiag;
    else
        Hest = complex(zeros(N,Nr,Nt,ns));
        for tx = 1:Nt
            pilotIndices = schedule.pilotIndicesByTx{tx};
            obs4 = reshape(observationsByTx{tx},N,Nr,1,numel(pilotIndices));
            prior4 = reshape(Hmean(:,:,tx,:),N,Nr,1,ns);
            if strcmp(method,'jakes-lmmse')
                pilotLSNoiseVariance = noiseVarianceFD/cfg.system.pilotPower;
                setup = temporal_channel_tracker_setup(pilotIndices,ns, ...
                    pilotLSNoiseVariance,params.maxDopplerHz,Tsym,diffuseVariance);
                [localEstimate,posteriorVariance] = temporal_channel_tracking_apply( ...
                    obs4,pilotIndices,ns,method,prior4,setup);
                posteriorMaps(:,:,tx,:) = repmat(reshape(posteriorVariance,1,1,1,ns),N,Nr,1,1);
            else
                localEstimate = temporal_channel_tracking_apply( ...
                    obs4,pilotIndices,ns,method,prior4,[]);
            end
            Hest(:,:,tx,:) = localEstimate;
        end
    end

    dataTruth = Hdiag(:,:,:,schedule.dataIndices);
    dataEstimate = Hest(:,:,:,schedule.dataIndices);
    channelError = dataEstimate-dataTruth;
    methods(im).name = method;
    methods(im).Hest = Hest;
    methods(im).MSEPerCoefficient = mean(abs(channelError(:)).^2);
    methods(im).NMSE = sum(abs(channelError(:)).^2)/max(sum(abs(dataTruth(:)).^2),eps);
    if strcmp(method,'jakes-lmmse')
        posteriorSelection = posteriorMaps(:,:,:,schedule.dataIndices);
        methods(im).PosteriorMSEPerCoefficient = mean(posteriorSelection(:));
        methods(im).PosteriorCalibrationRatio = methods(im).MSEPerCoefficient/ ...
            max(methods(im).PosteriorMSEPerCoefficient,eps);
    end

    bitErrors = 0;
    totalBits = numel(dataBits);
    for q = 1:schedule.numDataSymbols
        symbolIndex = schedule.dataIndices(q);
        Ysymbol = reshape(Yfd(:,:,symbolIndex),N,Nr);
        if strcmp(method,'perfect-csi-ici')
            Hfull = build_time_varying_ofdm_ici_matrix( ...
                tapSequence,N,cp,symbolIndex);
            [Xhat,~] = mimo_ici_equalize(Ysymbol,Hfull,noiseVarianceFD,'mmse');
        else
            Hsymbol = reshape(Hest(:,:,:,symbolIndex),N,Nr,Nt);
            [Xhat3,~] = mimo_equalize(reshape(Ysymbol,N,Nr,1), ...
                Hsymbol,noiseVarianceFD,'mmse');
            Xhat = Xhat3(:,:,1);
        end
        for tx = 1:Nt
            decodedBits = qam_demod(Xhat(:,tx),params.modulationOrder);
            bitErrors = bitErrors+nnz(decodedBits~=dataBits(:,tx,q));
        end
    end
    methods(im).BitErrors = bitErrors;
    methods(im).TotalBits = totalBits;
    methods(im).BER = bitErrors/totalBits;
    methods(im).FrameError = bitErrors>0;
end

frame.schedule = schedule;
frame.methods = methods;
frame.channelMeta = channelMeta;
frame.noiseVarianceFD = noiseVarianceFD;
frame.PilotContaminationVariance = pilotContaminationVariance;
frame.ICI_Ratio = iciRatio;
frame.PilotOverhead = schedule.pilotOverhead;
frame.NumPilotSymbols = schedule.numPilotSymbols;
frame.NumDataSymbols = schedule.numDataSymbols;
frame.NumBits = numel(dataBits);
frame.SpectralEfficiency = (schedule.numDataSymbols/ns)*Nt*log2(params.modulationOrder)*N/(N+cp);
frame.PilotIndicesByTx = schedule.pilotIndicesByTx;
frame.DataIndices = schedule.dataIndices;
end
