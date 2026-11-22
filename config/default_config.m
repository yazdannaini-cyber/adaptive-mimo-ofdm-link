function cfg = default_config()
%DEFAULT_CONFIG Configuration for the Adaptive MIMO-OFDM project.
% The project uses base MATLAB functions only (no Communications Toolbox).

cfg.system.nFFT = 64;
cfg.system.cpLength = 16;
cfg.system.nTx = 2;
cfg.system.nRx = 2;
cfg.system.nDataSymbols = 8;
cfg.system.pilotPower = 1.0;
cfg.system.subcarrierSpacingHz = 15e3; % 15 kHz reference numerology

cfg.channel.numTaps = 5;
cfg.channel.pdp = exp(-0.65 * (0:cfg.channel.numTaps-1));
cfg.channel.pdp = cfg.channel.pdp(:) / sum(cfg.channel.pdp);
cfg.channel.ricianK_dB = 6;

% Sample-varying Doppler channel model (V11).
cfg.channel.doppler.maxDopplerHz = 120;
cfg.channel.doppler.numOscillators = 64;
cfg.channel.doppler.losDopplerFraction = 1.0;

% Optional Kronecker spatial correlation for diffuse fading components.
% rho=0 gives independent antennas; larger rho increases coupling.
cfg.channel.spatialCorrelation.enabled = false;
cfg.channel.spatialCorrelation.txRho = 0;
cfg.channel.spatialCorrelation.rxRho = 0;

% Optional joint spatial-frequency LMMSE estimator (off by default).
cfg.estimation.spatialAwareLMMSE = false;

cfg.simulation.channelModels = {'awgn', 'rayleigh', 'rician'};
cfg.simulation.modulations = [4, 16, 64];
cfg.simulation.ebn0dB = 0:4:24;
cfg.simulation.numFrames = 80;
cfg.simulation.randomSeed = 20261121;

cfg.adaptive.snrDb = 0:4:28;          % Es/N0 per subcarrier, not Eb/N0.
cfg.adaptive.modulations = [4, 16, 64];
% Approximate uncoded square-QAM Es/N0 thresholds for a target BER near 1e-3.
% These are starting points; the resulting BER must be verified from simulation.
cfg.adaptive.thresholdDb = [9.8, 16.5, 22.5]; % retained as a selectable legacy policy
cfg.adaptive.selectionMethod = 'predicted-ber';
cfg.adaptive.targetBER = 1e-3;
cfg.adaptive.percentile = 10; % reporting metric and legacy threshold policy

% Frame-cluster bootstrap statistics (captures within-frame error dependence).
cfg.statistics.bootstrapSamples = 1000;
cfg.statistics.confidenceLevel = 0.95;
cfg.statistics.bootstrapSeed = 20261122;

cfg.output.directory = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
cfg.output.saveFigures = true;
cfg.output.saveMatFiles = true;

end
