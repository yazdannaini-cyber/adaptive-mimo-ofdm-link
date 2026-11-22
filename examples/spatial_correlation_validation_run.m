function results = spatial_correlation_validation_run()
%SPATIAL_CORRELATION_VALIDATION_RUN Stress-test adaptive MIMO under spatial correlation.
% From the project root:
%   addpath('examples');
%   results = spatial_correlation_validation_run;
%
% The current per-link LMMSE estimator does not exploit cross-antenna covariance.
% This is deliberately a robustness/sensitivity study of the V6 controller.
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
cfg = default_config();
cfg.adaptiveMimo.channelModels = {'rayleigh','rician'};
cfg.adaptiveMimo.snrDb = [16 24 32];
cfg.adaptiveMimo.modulations = [4 16 64];
cfg.adaptiveMimo.numFrames = 500;
cfg.adaptiveMimo.targetBER = 1e-3;
cfg.adaptiveMimo.bootstrapSamples = 1000;
cfg.adaptiveMimo.confidenceLevel = 0.95;
cfg.adaptiveMimo.bootstrapSeed = cfg.statistics.bootstrapSeed + 17000;

% Each row is [rhoTx rhoRx]. Separating TX-only and RX-only helps identify
% which side of the link causes the larger loss.
profiles = [0.0 0.0; 0.7 0.0; 0.0 0.7; 0.7 0.7; 0.9 0.9];
profileNames = ["iid";"tx-correlated";"rx-correlated";"joint-moderate";"joint-high"];
outDir = fullfile(projectRoot,'spatial_correlation_results');
if ~exist(outDir,'dir'), mkdir(outDir); end
rng(cfg.simulation.randomSeed+1700,'twister');
results = [];

for p = 1:size(profiles,1)
    fprintf('\n=== Correlation profile: %s (rhoTx=%.2f, rhoRx=%.2f) ===\n', ...
        char(profileNames(p)),profiles(p,1),profiles(p,2));
    cfg.channel.spatialCorrelation.enabled = true;
    cfg.channel.spatialCorrelation.txRho = profiles(p,1);
    cfg.channel.spatialCorrelation.rxRho = profiles(p,2);
    for c = 1:numel(cfg.adaptiveMimo.channelModels)
        channelName = cfg.adaptiveMimo.channelModels{c};
        current = run_adaptive_mimo_mode(cfg,channelName);
        current.TxCorrelationCoefficient = repmat(profiles(p,1),height(current),1);
        current.RxCorrelationCoefficient = repmat(profiles(p,2),height(current),1);
        current.CorrelationProfile = repmat(profileNames(p),height(current),1);
        if isempty(results), results = current; else, results = [results;current]; %#ok<AGROW>
        end
    end
end
writetable(results,fullfile(outDir,'spatial_correlation_validation.csv'));
save(fullfile(outDir,'spatial_correlation_validation.mat'),'cfg','profiles','profileNames','results');
plot_spatial_correlation_results(results,outDir);
fprintf('\nSPATIAL CORRELATION VALIDATION COMPLETED\n');
fprintf('Channels              : %s\n',strjoin(cfg.adaptiveMimo.channelModels,', '));
fprintf('Correlation profiles  : %d\n',size(profiles,1));
fprintf('SNR points/profile    : %d\n',numel(cfg.adaptiveMimo.snrDb));
fprintf('Frames/channel/SNR    : %d\n',cfg.adaptiveMimo.numFrames);
fprintf('Result rows           : %d\n',height(results));
fprintf('Output                : %s\n',outDir);
end
