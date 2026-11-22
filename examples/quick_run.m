function quick_run()
%QUICK_RUN Execute a small smoke test across all configured channel families.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
cfg = default_config();
cfg.simulation.channelModels = {'awgn','rayleigh','rician'};
cfg.simulation.ebn0dB = [0 8 16 24];
cfg.simulation.numFrames = 5;
cfg.statistics.bootstrapSamples = 300;
cfg.adaptive.snrDb = [0 8 16 24];
cfg.output.directory = fullfile(projectRoot,'quick_results');
run_project(cfg);
end
