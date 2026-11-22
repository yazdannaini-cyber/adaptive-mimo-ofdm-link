function results = mimo_mode_comparison_run()
%MIMO_MODE_COMPARISON_RUN Run the spatial-multiplexing vs Alamouti study.
% From the project root:
%   addpath('examples');
%   results = mimo_mode_comparison_run;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
cfg = default_config();
cfg.comparison.channelModels = {'awgn','rayleigh','rician'};
cfg.comparison.snrDb = 0:8:32;
cfg.comparison.modulations = [4 16 64];
cfg.comparison.numFrames = 200;
cfg.statistics.bootstrapSamples = 1000;
cfg.output.directory = fullfile(projectRoot,'comparison_results');
if ~exist(cfg.output.directory,'dir')
    mkdir(cfg.output.directory);
end
rng(cfg.simulation.randomSeed+500,'twister');
results = [];
for channelIndex = 1:numel(cfg.comparison.channelModels)
    channelName = cfg.comparison.channelModels{channelIndex};
    channelResults = run_mimo_mode_comparison(cfg,channelName);
    if isempty(results)
        results = channelResults;
    else
        results = [results;channelResults]; %#ok<AGROW>
    end
end
writetable(results,fullfile(cfg.output.directory,'mimo_mode_comparison_results.csv'));
save(fullfile(cfg.output.directory,'mimo_mode_comparison_results.mat'),'cfg','results');
plot_mimo_mode_comparison(results,cfg.output.directory);

fprintf('\nMIMO MODE COMPARISON COMPLETED\n');
fprintf('Channels              : %s\n',strjoin(cfg.comparison.channelModels,', '));
fprintf('Modulations           : %s\n',mat2str(cfg.comparison.modulations));
fprintf('SNR points            : %d\n',numel(cfg.comparison.snrDb));
fprintf('Frames/configuration  : %d\n',cfg.comparison.numFrames);
fprintf('Result rows           : %d\n',height(results));
fprintf('Output                 : %s\n',cfg.output.directory);
end
