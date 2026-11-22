function results = spatial_aware_lmmse_comparison_run()
%SPATIAL_AWARE_LMMSE_COMPARISON_RUN Compare per-link and joint LMMSE fairly.
% Each pair uses deterministic frame-level seeds, so channel and pilot-noise
% realizations are exactly matched between estimators for each profile/channel.
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
cfg.adaptiveMimo.bootstrapSeed = cfg.statistics.bootstrapSeed + 28000;

profiles = [0.0 0.0; 0.7 0.0; 0.0 0.7; 0.7 0.7; 0.9 0.9];
profileNames = ["iid";"tx-correlated";"rx-correlated";"joint-moderate";"joint-high"];
estimatorFlags = [false true];
estimatorNames = ["per-link-lmmse";"joint-spatial-frequency-lmmse"];
outDir = fullfile(projectRoot,'spatial_aware_lmmse_results');
if ~exist(outDir,'dir'), mkdir(outDir); end
results = [];

for p = 1:size(profiles,1)
    cfg.channel.spatialCorrelation.enabled = true;
    cfg.channel.spatialCorrelation.txRho = profiles(p,1);
    cfg.channel.spatialCorrelation.rxRho = profiles(p,2);
    fprintf('\n============================================================\n');
    fprintf('PROFILE %s | rhoTx=%.2f | rhoRx=%.2f\n', ...
        char(profileNames(p)),profiles(p,1),profiles(p,2));
    for c = 1:numel(cfg.adaptiveMimo.channelModels)
        channelName = cfg.adaptiveMimo.channelModels{c};
        for e = 1:numel(estimatorFlags)
            cfg.estimation.spatialAwareLMMSE = estimatorFlags(e);
            fprintf('\nEstimator: %s | Channel: %s\n', ...
                char(estimatorNames(e)),upper(char(channelName)));
            % Paired channel/pilot randomness: the intended difference is the estimator.
            pairSeed = cfg.simulation.randomSeed + 40000 + 1000*p + 100*c;
            cfg.adaptiveMimo.frameSeedBase = pairSeed;
            rng(pairSeed,'twister');
            current = run_adaptive_mimo_mode(cfg,channelName);
            current.Estimator = repmat(estimatorNames(e),height(current),1);
            current.CorrelationProfile = repmat(profileNames(p),height(current),1);
            current.TxCorrelationCoefficient = repmat(profiles(p,1),height(current),1);
            current.RxCorrelationCoefficient = repmat(profiles(p,2),height(current),1);
            if isempty(results), results = current; else, results = [results;current]; %#ok<AGROW>
            end
        end
    end
end
writetable(results,fullfile(outDir,'spatial_aware_lmmse_comparison.csv'));
save(fullfile(outDir,'spatial_aware_lmmse_comparison.mat'), ...
    'cfg','profiles','profileNames','estimatorNames','results');
plot_spatial_aware_lmmse_comparison(results,outDir);
fprintf('\nSPATIAL-AWARE LMMSE COMPARISON COMPLETED\n');
fprintf('Channels             : %s\n',strjoin(cfg.adaptiveMimo.channelModels,', '));
fprintf('Correlation profiles : %d\n',size(profiles,1));
fprintf('Estimators           : %s\n',strjoin(cellstr(estimatorNames),', '));
fprintf('SNR points/profile   : %d\n',numel(cfg.adaptiveMimo.snrDb));
fprintf('Frames/run           : %d\n',cfg.adaptiveMimo.numFrames);
fprintf('Result rows          : %d\n',height(results));
fprintf('Output               : %s\n',outDir);
end
