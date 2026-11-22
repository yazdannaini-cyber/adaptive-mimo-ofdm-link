function validate_config(cfg)
%VALIDATE_CONFIG Fail early when simulation settings are inconsistent.
requiredTopLevel = {'system','channel','simulation','adaptive','output','statistics'};
for k = 1:numel(requiredTopLevel)
    if ~isfield(cfg,requiredTopLevel{k})
        error('validate_config:MissingField','Missing cfg.%s.',requiredTopLevel{k});
    end
end

if cfg.system.nFFT < 16 || cfg.system.nFFT ~= floor(cfg.system.nFFT)
    error('validate_config:FFT','nFFT must be an integer >= 16.');
end
if cfg.system.nTx < 1 || cfg.system.nRx < cfg.system.nTx || ...
        cfg.system.nTx ~= floor(cfg.system.nTx) || cfg.system.nRx ~= floor(cfg.system.nRx)
    error('validate_config:MIMO','This project requires nRx >= nTx >= 1.');
end
if cfg.channel.numTaps < 1 || cfg.channel.numTaps ~= floor(cfg.channel.numTaps) || ...
        cfg.system.cpLength ~= floor(cfg.system.cpLength) || ...
        cfg.system.cpLength < cfg.channel.numTaps-1 || cfg.system.cpLength >= cfg.system.nFFT
    error('validate_config:CP','CP must be >= channel memory and < nFFT.');
end
if cfg.system.nDataSymbols < 1 || cfg.simulation.numFrames < 1 || ...
        cfg.system.nDataSymbols ~= floor(cfg.system.nDataSymbols) || ...
        cfg.simulation.numFrames ~= floor(cfg.simulation.numFrames)
    error('validate_config:Frames','Number of data symbols and frames must be positive.');
end
if cfg.system.pilotPower <= 0
    error('validate_config:PilotPower','Pilot power must be positive.');
end
if numel(cfg.channel.pdp) ~= cfg.channel.numTaps || any(cfg.channel.pdp < 0) || ...
        abs(sum(cfg.channel.pdp)-1) > 1e-8
    error('validate_config:PDP','PDP must be nonnegative, have numTaps elements, and sum to one.');
end
for M = cfg.simulation.modulations
    if M < 4 || abs(log2(M)-round(log2(M))) > eps || ...
            abs(log2(sqrt(M))-round(log2(sqrt(M)))) > eps
        error('validate_config:QAM','Simulation modulation orders must be square powers of four.');
    end
end
if isempty(cfg.simulation.channelModels) || isempty(cfg.simulation.modulations) || ...
        isempty(cfg.simulation.ebn0dB) || isempty(cfg.adaptive.snrDb)
    error('validate_config:EmptySweep','Channel models, modulation orders, and SNR sweeps must be nonempty.');
end
if numel(cfg.adaptive.modulations) ~= numel(cfg.adaptive.thresholdDb) || ...
        any(diff(cfg.adaptive.thresholdDb) <= 0) || any(diff(cfg.adaptive.modulations) <= 0)
    error('validate_config:Adaptive','Adaptive modulation orders and thresholds must be strictly increasing and paired.');
end
if numel(cfg.adaptive.modulations) ~= 3 || ~isequal(cfg.adaptive.modulations,[4 16 64])
    error('validate_config:AdaptiveModes','Current adaptive result schema expects [4 16 64] QAM modes.');
end
if cfg.adaptive.percentile <= 0 || cfg.adaptive.percentile > 100
    error('validate_config:Percentile','Adaptive percentile must be in (0,100].');
end
selectionMethod = lower(char(cfg.adaptive.selectionMethod));
if ~strcmp(selectionMethod,'predicted-ber') && ~strcmp(selectionMethod,'sinr-threshold')
    error('validate_config:SelectionMethod', ...
        'Adaptive selection method must be predicted-ber or sinr-threshold.');
end
if ~isscalar(cfg.adaptive.targetBER) || (isnan(cfg.adaptive.targetBER) || isinf(cfg.adaptive.targetBER)) || ...
        cfg.adaptive.targetBER <= 0 || cfg.adaptive.targetBER >= 0.5
    error('validate_config:TargetBER','Adaptive target BER must be in (0,0.5).');
end
if ~isscalar(cfg.statistics.bootstrapSamples) || cfg.statistics.bootstrapSamples < 200 || ...
        cfg.statistics.bootstrapSamples ~= floor(cfg.statistics.bootstrapSamples) || ...
        ~isscalar(cfg.statistics.confidenceLevel) || cfg.statistics.confidenceLevel <= 0 || ...
        cfg.statistics.confidenceLevel >= 1 || ~isscalar(cfg.statistics.bootstrapSeed) || ...
        cfg.statistics.bootstrapSeed < 0 || cfg.statistics.bootstrapSeed ~= floor(cfg.statistics.bootstrapSeed)
    error('validate_config:Statistics','Invalid bootstrap settings. Use >=200 samples, confidence in (0,1), and a nonnegative integer seed.');
end

if isfield(cfg.system,'subcarrierSpacingHz') && ...
        (~isscalar(cfg.system.subcarrierSpacingHz) || ~isreal(cfg.system.subcarrierSpacingHz) || ...
        ~isfinite(cfg.system.subcarrierSpacingHz) || cfg.system.subcarrierSpacingHz <= 0)
    error('validate_config:SubcarrierSpacing','subcarrierSpacingHz must be finite and positive.');
end
if isfield(cfg.channel,'doppler')
    d = cfg.channel.doppler;
    if ~all(isfield(d,{'maxDopplerHz','numOscillators','losDopplerFraction'})) || ...
            ~isscalar(d.maxDopplerHz) || ~isreal(d.maxDopplerHz) || ...
            ~isfinite(d.maxDopplerHz) || d.maxDopplerHz < 0 || ...
            ~isscalar(d.numOscillators) || ~isreal(d.numOscillators) || ...
            ~isfinite(d.numOscillators) || d.numOscillators < 8 || ...
            d.numOscillators ~= floor(d.numOscillators) || ...
            ~isscalar(d.losDopplerFraction) || ~isreal(d.losDopplerFraction) || ...
            ~isfinite(d.losDopplerFraction) || abs(d.losDopplerFraction) > 1
        error('validate_config:Doppler','Invalid Doppler channel configuration.');
    end
end

if isfield(cfg,'estimation') && isfield(cfg.estimation,'spatialAwareLMMSE')
    flag = cfg.estimation.spatialAwareLMMSE;
    validFlag = isscalar(flag) && (islogical(flag) || ...
        (isnumeric(flag) && isreal(flag) && isfinite(flag) && any(flag == [0 1])));
    if ~validFlag
        error('validate_config:SpatialAwareLMMSE', ...
            'cfg.estimation.spatialAwareLMMSE must be a logical scalar.');
    end
end

end
