function results = adaptive_mimo_mode_run()
%ADAPTIVE_MIMO_MODE_RUN Run joint MIMO mode + QAM adaptation validation.
% From the project root:
%   addpath('examples');
%   results = adaptive_mimo_mode_run;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
cfg = default_config();
cfg.adaptiveMimo.channelModels = {'rayleigh','rician'};
cfg.adaptiveMimo.snrDb = [12 16 20 24 28 32];
cfg.adaptiveMimo.modulations = [4 16 64];
cfg.adaptiveMimo.numFrames = 1000;
cfg.adaptiveMimo.targetBER = 1e-3;
cfg.adaptiveMimo.bootstrapSamples = 2000;
cfg.adaptiveMimo.confidenceLevel = 0.95;
cfg.adaptiveMimo.bootstrapSeed = cfg.statistics.bootstrapSeed + 9000;
cfg.output.directory = fullfile(projectRoot,'adaptive_mimo_results');
if ~exist(cfg.output.directory,'dir'), mkdir(cfg.output.directory); end
rng(cfg.simulation.randomSeed+700,'twister');
results = [];
for c = 1:numel(cfg.adaptiveMimo.channelModels)
    current = run_adaptive_mimo_mode(cfg,cfg.adaptiveMimo.channelModels{c});
    if isempty(results), results = current; else, results = [results;current]; %#ok<AGROW>
    end
end
writetable(results,fullfile(cfg.output.directory,'adaptive_mimo_mode_results.csv'));
save(fullfile(cfg.output.directory,'adaptive_mimo_mode_results.mat'),'cfg','results');
plot_adaptive_mimo_mode_results(results,cfg.output.directory);
fprintf('\nADAPTIVE MIMO MODE STUDY COMPLETED\n');
fprintf('Channels             : %s\n',strjoin(cfg.adaptiveMimo.channelModels,', '));
fprintf('SNR points           : %d\n',numel(cfg.adaptiveMimo.snrDb));
fprintf('Frames/channel/SNR   : %d\n',cfg.adaptiveMimo.numFrames);
fprintf('Result rows          : %d\n',height(results));
fprintf('Output               : %s\n',cfg.output.directory);
end
