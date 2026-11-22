function [fixedResults, adaptiveResults] = run_project(cfg)
%RUN_PROJECT Execute the complete MIMO-OFDM simulation study.
% Run this file from MATLAB: >> run_project

projectRoot = fileparts(mfilename('fullpath'));
addpath(genpath(projectRoot));
if nargin < 1 || isempty(cfg)
    cfg = default_config();
end
validate_config(cfg);
rng(cfg.simulation.randomSeed,'twister');
outputDir = cfg.output.directory;
if ~exist(outputDir,'dir'), mkdir(outputDir); end

fixedResults = [];
adaptiveResults = [];
for c = 1:numel(cfg.simulation.channelModels)
    channelName = cfg.simulation.channelModels{c};
    fixedPart = run_fixed_modulation(cfg,channelName);
    adaptivePart = run_adaptive_modulation(cfg,channelName);
    if c == 1
        fixedResults = fixedPart;
        adaptiveResults = adaptivePart;
    else
        fixedResults = [fixedResults; fixedPart]; %#ok<AGROW>
        adaptiveResults = [adaptiveResults; adaptivePart]; %#ok<AGROW>
    end
end

writetable(fixedResults,fullfile(outputDir,'fixed_modulation_results.csv'));
writetable(adaptiveResults,fullfile(outputDir,'adaptive_modulation_results.csv'));
if cfg.output.saveMatFiles
    save(fullfile(outputDir,'complete_results.mat'),'cfg','fixedResults','adaptiveResults','-v7.3');
end
plot_results(fixedResults,adaptiveResults,outputDir,cfg.output.saveFigures);

fprintf('\n============================================================\n');
fprintf('MIMO-OFDM SIMULATION COMPLETED\n');
fprintf('Fixed-modulation rows : %d\n',height(fixedResults));
fprintf('Adaptive-modulation rows: %d\n',height(adaptiveResults));
fprintf('Results directory     : %s\n',outputDir);
fprintf('============================================================\n');
end
