function results = adaptive_validation_run()
%ADAPTIVE_VALIDATION_RUN Longer Monte-Carlo study focused on adaptive link reliability.
% Run from the project root:
%   addpath('examples');
%   results = adaptive_validation_run;
%
% This intentionally runs adaptive simulations only (not the much larger fixed
% modulation sweep). Increase numFrames further for final publication-quality
% results if the BER confidence intervals remain too wide.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
cfg = default_config();
cfg.simulation.channelModels = {'rayleigh','rician'};
cfg.simulation.numFrames = 1000;
cfg.adaptive.snrDb = [12 16 20 24 28 32];
cfg.statistics.bootstrapSamples = 2000;
cfg.output.directory = fullfile(projectRoot,'validation_results');
if ~exist(cfg.output.directory,'dir')
    mkdir(cfg.output.directory);
end
validate_config(cfg);
rng(cfg.simulation.randomSeed,'twister');
results = [];
for c = 1:numel(cfg.simulation.channelModels)
    channelName = cfg.simulation.channelModels{c};
    part = run_adaptive_modulation(cfg,channelName);
    if isempty(results)
        results = part;
    else
        results = [results;part]; %#ok<AGROW>
    end
end
writetable(results,fullfile(cfg.output.directory,'adaptive_validation_results.csv'));
save(fullfile(cfg.output.directory,'adaptive_validation_results.mat'),'cfg','results');

fprintf('\nADAPTIVE VALIDATION STUDY COMPLETED\n');
fprintf('Frames per SNR/channel : %d\n',cfg.simulation.numFrames);
fprintf('Bootstrap resamples    : %d\n',cfg.statistics.bootstrapSamples);
fprintf('Output                 : %s\n',cfg.output.directory);
end
