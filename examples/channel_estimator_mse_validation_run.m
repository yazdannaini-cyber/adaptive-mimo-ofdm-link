function results = channel_estimator_mse_validation_run()
%CHANNEL_ESTIMATOR_MSE_VALIDATION_RUN Direct paired NMSE study of LMMSE estimators.
% Measures channel-estimation error directly, independently of adaptive
% modulation or MIMO mode selection. Outputs paired confidence intervals and
% posterior-MSE calibration metrics for both Rayleigh and Rician channels.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
cfg = default_config();
cfg.channel.spatialCorrelation.enabled = true;
profiles = [0.0 0.0;0.7 0.0;0.0 0.7;0.7 0.7;0.9 0.9];
profileNames = ["iid";"tx-correlated";"rx-correlated";"joint-moderate";"joint-high"];
channelNames = {'rayleigh','rician'};
pilotSNRdB = [0 8 16 24 32];
numFrames = 500;
numBootstrap = 1000;
confidenceLevel = 0.95;
baseSeed = cfg.simulation.randomSeed+910000;

outDir = fullfile(projectRoot,'estimator_mse_results');
if ~exist(outDir,'dir'), mkdir(outDir); end
results = [];
for p = 1:size(profiles,1)
    cfg.channel.spatialCorrelation.txRho = profiles(p,1);
    cfg.channel.spatialCorrelation.rxRho = profiles(p,2);
    fprintf('\n============================================================\n');
    fprintf('PROFILE %s | rhoTx=%.2f | rhoRx=%.2f\n', ...
        char(profileNames(p)),profiles(p,1),profiles(p,2));
    for c = 1:numel(channelNames)
        for s = 1:numel(pilotSNRdB)
            seed = baseSeed+p*100000+c*10000+s*1000;
            row = channel_estimator_mse_pair(cfg,channelNames{c},pilotSNRdB(s), ...
                numFrames,numBootstrap,confidenceLevel,seed);
            row.CorrelationProfile = repmat(profileNames(p),height(row),1);
            row = movevars(row,'CorrelationProfile','After','Channel');
            if isempty(results), results = row; else, results = [results;row]; %#ok<AGROW>
            end
            fprintf(['  %s | pilot SNR=%2d dB | per-link NMSE=%.4g (%.2f dB) | ' ...
                'joint NMSE=%.4g (%.2f dB) | improvement=%+.1f%% | ' ...
                'paired MSE diff CI=[%.3g, %.3g]\n'], ...
                upper(channelNames{c}),pilotSNRdB(s),row.PerLink_NMSE,row.PerLink_NMSE_dB, ...
                row.Joint_NMSE,row.Joint_NMSE_dB,row.Joint_NMSE_Improvement_Percent, ...
                row.Paired_MSE_Improvement_CI95_Lower,row.Paired_MSE_Improvement_CI95_Upper);
        end
    end
end
writetable(results,fullfile(outDir,'channel_estimator_mse_validation.csv'));
save(fullfile(outDir,'channel_estimator_mse_validation.mat'), ...
    'cfg','profiles','profileNames','channelNames','pilotSNRdB','numFrames', ...
    'numBootstrap','confidenceLevel','results');
plot_channel_estimator_mse_validation(results,outDir);
fprintf('\nCHANNEL ESTIMATOR MSE VALIDATION COMPLETED\n');
fprintf('Channels              : %s\n',strjoin(channelNames,', '));
fprintf('Correlation profiles  : %d\n',size(profiles,1));
fprintf('Pilot SNR points      : %d\n',numel(pilotSNRdB));
fprintf('Estimators            : per-link LMMSE, joint spatial-frequency LMMSE\n');
fprintf('Frames/condition      : %d\n',numFrames);
fprintf('Bootstrap resamples   : %d\n',numBootstrap);
fprintf('Result rows           : %d\n',height(results));
fprintf('Output                : %s\n',outDir);
end
