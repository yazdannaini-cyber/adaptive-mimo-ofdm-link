function [tapSequence, components, meta] = generate_time_varying_mimo_channel( ...
    cfg, modelName, nSamples, sampleRateHz, seed)
%GENERATE_TIME_VARYING_MIMO_CHANNEL Generate sample-varying MIMO FIR taps.
% tapSequence dimensions: [nSamples x nTaps x nRx x nTx].
%
% The diffuse fading process uses a finite sum-of-sinusoids approximation
% to the classical isotropic-scattering/Jakes Doppler spectrum:
%   h(t) = 1/sqrt(Nosc) sum_m exp(j(2*pi*fD*cos(theta_m)*t + phi_m)).
% Independent diffuse links are then shaped by the configured Kronecker
% transmit/receive spatial correlation factors. Rician channels add a
% deterministic first-tap LOS component with a configurable Doppler shift.
%
% This is a discrete-time, sample-varying FIR channel. Because coefficients
% may change within an OFDM symbol, the channel can create inter-carrier
% interference (ICI); it is not a block-static per-symbol approximation.
%
% Inputs:
%   cfg          project configuration
%   modelName    'rayleigh', 'rician', or 'awgn'
%   nSamples     number of channel time samples (positive integer)
%   sampleRateHz channel sample rate in Hz
%   seed         optional nonnegative integer for reproducible generation
%
% Outputs:
%   tapSequence  time-varying channel taps
%   components   struct with separately scaled diffuse and LOS components
%   meta         Doppler and sampling metadata

if nargin < 4 || nargin > 5
    error('generate_time_varying_mimo_channel:InputCount', ...
        'Provide cfg, modelName, nSamples, sampleRateHz, and optional seed.');
end
if ~isstruct(cfg) || ~isfield(cfg,'system') || ~isfield(cfg,'channel') || ...
        ~isfield(cfg.channel,'pdp') || ~isfield(cfg.channel,'numTaps')
    error('generate_time_varying_mimo_channel:InvalidConfig', ...
        'cfg must contain system and channel configuration fields.');
end
if ~isscalar(nSamples) || ~isreal(nSamples) || ~isfinite(nSamples) || ...
        nSamples < 1 || nSamples ~= floor(nSamples)
    error('generate_time_varying_mimo_channel:InvalidSampleCount', ...
        'nSamples must be a positive integer.');
end
if ~isscalar(sampleRateHz) || ~isreal(sampleRateHz) || ...
        ~isfinite(sampleRateHz) || sampleRateHz <= 0
    error('generate_time_varying_mimo_channel:InvalidSampleRate', ...
        'sampleRateHz must be finite and positive.');
end
if nargin >= 5 && (~isscalar(seed) || ~isreal(seed) || ~isfinite(seed) || ...
        seed < 0 || seed ~= floor(seed))
    error('generate_time_varying_mimo_channel:InvalidSeed', ...
        'seed must be a nonnegative integer.');
end

modelName = lower(char(modelName));
if ~any(strcmp(modelName,{'rayleigh','rician','awgn'}))
    error('generate_time_varying_mimo_channel:UnknownModel', ...
        'modelName must be rayleigh, rician, or awgn.');
end

nRx = cfg.system.nRx;
nTx = cfg.system.nTx;
nTaps = cfg.channel.numTaps;
pdp = cfg.channel.pdp(:);
if numel(pdp) ~= nTaps || any(~isfinite(pdp)) || any(pdp < 0) || sum(pdp) <= 0
    error('generate_time_varying_mimo_channel:InvalidPDP', ...
        'PDP must contain numTaps finite nonnegative values with positive sum.');
end
pdp = pdp/sum(pdp);
if cfg.system.cpLength < nTaps-1
    error('generate_time_varying_mimo_channel:CPTooShort', ...
        'Cyclic prefix must be at least nTaps-1.');
end

if isfield(cfg.channel,'doppler')
    dopplerCfg = cfg.channel.doppler;
else
    dopplerCfg.maxDopplerHz = 0;
    dopplerCfg.numOscillators = 64;
    dopplerCfg.losDopplerFraction = 1;
end
requiredDopplerFields = {'maxDopplerHz','numOscillators','losDopplerFraction'};
if ~all(isfield(dopplerCfg,requiredDopplerFields))
    error('generate_time_varying_mimo_channel:MissingDopplerConfig', ...
        'Doppler config must define maxDopplerHz, numOscillators, and losDopplerFraction.');
end
fd = dopplerCfg.maxDopplerHz;
numOsc = dopplerCfg.numOscillators;
losFraction = dopplerCfg.losDopplerFraction;
if ~isscalar(fd) || ~isreal(fd) || ~isfinite(fd) || fd < 0 || fd >= sampleRateHz/2
    error('generate_time_varying_mimo_channel:InvalidDoppler', ...
        'maxDopplerHz must be in [0, sampleRateHz/2).');
end
if ~isscalar(numOsc) || ~isreal(numOsc) || ~isfinite(numOsc) || ...
        numOsc < 8 || numOsc ~= floor(numOsc)
    error('generate_time_varying_mimo_channel:InvalidOscillatorCount', ...
        'numOscillators must be an integer >= 8.');
end
if ~isscalar(losFraction) || ~isreal(losFraction) || ~isfinite(losFraction) || ...
        abs(losFraction) > 1
    error('generate_time_varying_mimo_channel:InvalidLOSDoppler', ...
        'losDopplerFraction must be a finite scalar in [-1,1].');
end

% Preserve the caller's global random stream when a local seed is supplied.
if nargin >= 5
    previousRng = rng;
    cleanupRng = onCleanup(@() rng(previousRng)); %#ok<NASGU>
    rng(seed,'twister');
end

if strcmp(modelName,'awgn')
    if nRx < nTx
        error('generate_time_varying_mimo_channel:AWGNTopology', ...
            'AWGN reference channel requires nRx >= nTx.');
    end
    tapSequence = complex(zeros(nSamples,nTaps,nRx,nTx));
    for rx = 1:nRx
        for tx = 1:nTx
            tapSequence(:,1,rx,tx) = double(rx == tx);
        end
    end
    components.diffuse = complex(zeros(size(tapSequence)));
    components.los = tapSequence;
    meta = makeMeta(modelName,fd,sampleRateHz,numOsc,nSamples,0);
    return;
end

[Lrx,Ltx,Rrx,Rtx] = spatial_correlation_factors(cfg);
time = (0:nSamples-1).'/sampleRateHz;
diffuseRaw = complex(zeros(nSamples,nTaps,nRx,nTx));
for ell = 1:nTaps
    for rx = 1:nRx
        for tx = 1:nTx
            theta = 2*pi*rand(numOsc,1);
            phase = 2*pi*rand(numOsc,1);
            toneHz = fd*cos(theta);
            phaseMatrix = 2*pi*(time*toneHz.') + phase.';
            process = sum(exp(1i*phaseMatrix),2)/sqrt(numOsc);
            diffuseRaw(:,ell,rx,tx) = process;
        end
    end
end

% Apply the separable Kronecker spatial model at every time sample and tap.
diffuseCorrelated = complex(zeros(size(diffuseRaw)));
for sample = 1:nSamples
    for ell = 1:nTaps
        Hraw = reshape(diffuseRaw(sample,ell,:,:),nRx,nTx);
        Hcorr = Lrx*Hraw*Ltx';
        diffuseCorrelated(sample,ell,:,:) = reshape(Hcorr,1,1,nRx,nTx);
    end
end

components.diffuse = complex(zeros(size(diffuseCorrelated)));
components.los = complex(zeros(size(diffuseCorrelated)));
if strcmp(modelName,'rayleigh')
    for ell = 1:nTaps
        components.diffuse(:,ell,:,:) = sqrt(pdp(ell))*diffuseCorrelated(:,ell,:,:);
    end
    K = NaN;
    losDopplerHz = 0;
else
    K = 10^(cfg.channel.ricianK_dB/10);
    diffuseScale = sqrt(1/(K+1));
    losScale = sqrt(K/(K+1));
    for ell = 1:nTaps
        components.diffuse(:,ell,:,:) = sqrt(pdp(ell))*diffuseScale*diffuseCorrelated(:,ell,:,:);
    end
    losArray = complex(zeros(nRx,nTx));
    for rx = 1:nRx
        for tx = 1:nTx
            losArray(rx,tx) = exp(-1i*2*pi*(rx-1)*(tx-1)/max(nRx,nTx));
        end
    end
    losDopplerHz = fd*losFraction;
    losTime = exp(1i*2*pi*losDopplerHz*time);
    for rx = 1:nRx
        for tx = 1:nTx
            components.los(:,1,rx,tx) = sqrt(pdp(1))*losScale* ...
                losArray(rx,tx)*losTime;
        end
    end
end

tapSequence = components.diffuse + components.los;
meta = makeMeta(modelName,fd,sampleRateHz,numOsc,nSamples,losDopplerHz);
meta.txCorrelation = Rtx;
meta.rxCorrelation = Rrx;
meta.pdp = pdp;
meta.ricianK_dB = cfg.channel.ricianK_dB;
end

function meta = makeMeta(modelName,fd,sampleRateHz,numOsc,nSamples,losDopplerHz)
meta.modelName = modelName;
meta.maxDopplerHz = fd;
meta.sampleRateHz = sampleRateHz;
meta.numOscillators = numOsc;
meta.numSamples = nSamples;
meta.losDopplerHz = losDopplerHz;
if fd > 0
    meta.coherenceTimeApprox_s = 0.423/fd;
else
    meta.coherenceTimeApprox_s = Inf;
end
end
