function results = pilot_mapped_adaptive_validation_run()
%PILOT_MAPPED_ADAPTIVE_VALIDATION_RUN V15 frame-level link-adaptation study.
% Runs an actual pilot-mapped, sample-varying 2x2 MIMO-OFDM waveform. The
% transmitter chooses between Spatial Multiplexing and Alamouti STBC, and
% between QPSK and 16-QAM, using only the first orthogonal pilot occasion.
% Later pilot observations feed the causal Jakes-LMMSE data-channel tracker.
%
% Output: pilot_mapped_results_v15/adaptive_link_results.csv and MAT/PNG files.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
cfg = default_config();
cfg.channel.doppler.numOscillators = 32;

channelNames = {'rayleigh','rician'};
profileNames = {'iid','joint-high'};
profileRhos = [0 0; 0.9 0.9];
dopplerValuesHz = [0 120 300];
pilotSNRValuesdB = [8 14 20];
pilotSpacing = 8;
nSymbols = 24;
modulations = [4 16];
targetBER = 1e-2;
numFrames = 8;       % Screening run; raise to 100+ for publication-quality confidence.
numBootstrap = 1000;
confidenceLevel = 0.95;
strategies = {'adaptive','fixed-sm-16qam','fixed-alamouti-qpsk'};

outDir = fullfile(projectRoot,'pilot_mapped_results_v15');
if ~exist(outDir,'dir'), mkdir(outDir); end
results = [];
frameLog = [];
conditionCount = numel(channelNames)*numel(profileNames)* ...
    numel(dopplerValuesHz)*numel(pilotSNRValuesdB);
conditionNumber = 0;
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>

for c = 1:numel(channelNames)
    for p = 1:numel(profileNames)
        for d = 1:numel(dopplerValuesHz)
            for is = 1:numel(pilotSNRValuesdB)
                conditionNumber = conditionNumber+1;
                channelName = channelNames{c};
                profileName = profileNames{p};
                fd = dopplerValuesHz(d);
                snrDb = pilotSNRValuesdB(is);
                rhoTx = profileRhos(p,1);
                rhoRx = profileRhos(p,2);
                fprintf('\n[%d/%d] %s | %s | fD=%g Hz | pilot/data SNR=%g dB\n', ...
                    conditionNumber,conditionCount,upper(channelName),profileName,fd,snrDb);

                frameErrors = zeros(numFrames,numel(strategies));
                frameBits = zeros(numFrames,numel(strategies));
                frameBER = zeros(numFrames,numel(strategies));
                frameFER = false(numFrames,numel(strategies));
                frameSE = zeros(numFrames,numel(strategies));
                frameGoodput = zeros(numFrames,numel(strategies));
                frameNMSE = zeros(numFrames,numel(strategies));
                frameICI = zeros(numFrames,numel(strategies));
                framePredictedBER = zeros(numFrames,numel(strategies));
                frameFeasible = false(numFrames,numel(strategies));
                modeIndex = zeros(numFrames,1);
                modulationOrder = zeros(numFrames,1);

                for frameIndex = 1:numFrames
                    seedBase = cfg.simulation.randomSeed + 1e7*c + 1e6*p + ...
                        1e5*d + 1e4*is + 100*frameIndex;
                    params.nSymbols = nSymbols;
                    params.pilotSpacing = pilotSpacing;
                    params.pilotSNRdB = snrDb;
                    params.maxDopplerHz = fd;
                    params.txRho = rhoTx;
                    params.rxRho = rhoRx;
                    params.targetBER = targetBER;
                    params.modulations = modulations;
                    params.estimator = 'jakes-lmmse';
                    params.dataSeed = seedBase+1;
                    params.pilotSeed = seedBase+2;
                    params.channelSeed = seedBase+3;
                    params.noiseSeed = seedBase+4;

                    frame = simulate_pilot_mapped_adaptive_frame(cfg,channelName,params);
                    candidates = {frame.adaptive,frame.fixedSM16QAM,frame.fixedAlamoutiQPSK};
                    for strategyIndex = 1:numel(strategies)
                        candidate = candidates{strategyIndex};
                        frameErrors(frameIndex,strategyIndex) = candidate.BitErrors;
                        frameBits(frameIndex,strategyIndex) = candidate.TotalBits;
                        frameBER(frameIndex,strategyIndex) = candidate.BER;
                        frameFER(frameIndex,strategyIndex) = candidate.FrameError;
                        frameSE(frameIndex,strategyIndex) = candidate.NetSpectralEfficiency_bpsHz;
                        frameGoodput(frameIndex,strategyIndex) = candidate.UncodedGoodput_bpsHz;
                        frameNMSE(frameIndex,strategyIndex) = candidate.NMSE;
                        frameICI(frameIndex,strategyIndex) = candidate.ICI_Ratio;
                        framePredictedBER(frameIndex,strategyIndex) = candidate.PredictedBER;
                        frameFeasible(frameIndex,strategyIndex) = candidate.TargetFeasible;
                    end
                    modeIndex(frameIndex) = double(strcmp(frame.adaptive.Mode,'alamouti-stbc'))+1;
                    modulationOrder(frameIndex) = frame.adaptive.ModulationOrder;
                    logRow = table(string(channelName),string(profileName),fd,snrDb, ...
                        frameIndex,string(frame.adaptive.Mode),frame.adaptive.ModulationOrder, ...
                        frame.adaptive.PredictedBER,frame.adaptive.TargetFeasible, ...
                        frame.adaptive.BER,frame.adaptive.NetSpectralEfficiency_bpsHz, ...
                        frame.fixedSM16QAM.BER,frame.fixedAlamoutiQPSK.BER, ...
                        frame.adaptive.NMSE,frame.adaptive.ICI_Ratio, ...
                        'VariableNames',{'Channel','SpatialProfile','MaxDoppler_Hz','PilotSNR_dB', ...
                        'FrameIndex','SelectedMode','SelectedModulationOrder','PredictedBER', ...
                        'TargetFeasible','AdaptiveBER','AdaptiveNetSE_bpsHz', ...
                        'FixedSM16QAM_BER','FixedAlamoutiQPSK_BER','AdaptiveNMSE','AdaptiveICI_Ratio'});
                    if isempty(frameLog), frameLog = logRow; else, frameLog = [frameLog;logRow]; end %#ok<AGROW>
                end

                for strategyIndex = 1:numel(strategies)
                    totalErrors = sum(frameErrors(:,strategyIndex));
                    totalBits = sum(frameBits(:,strategyIndex));
                    totalFrameErrors = nnz(frameFER(:,strategyIndex));
                    metrics = bit_error_metrics(totalErrors,totalBits,totalFrameErrors,numFrames);
                    bootstrapCI = cluster_bootstrap_ber_ci(frameErrors(:,strategyIndex), ...
                        frameBits(:,strategyIndex),numBootstrap,confidenceLevel, ...
                        cfg.statistics.bootstrapSeed+conditionNumber*100+strategyIndex);
                    if strategyIndex == 1
                        smRatio = mean(modeIndex==1);
                        alamoutiRatio = mean(modeIndex==2);
                        qpskRatio = mean(modulationOrder==4);
                        qam16Ratio = mean(modulationOrder==16);
                    elseif strategyIndex == 2
                        smRatio = 1; alamoutiRatio = 0; qpskRatio = 0; qam16Ratio = 1;
                    else
                        smRatio = 0; alamoutiRatio = 1; qpskRatio = 1; qam16Ratio = 0;
                    end
                    row = table(string(channelName),string(profileName),rhoTx,rhoRx,fd,snrDb, ...
                        pilotSpacing,targetBER,string(strategies{strategyIndex}), ...
                        mean(framePredictedBER(:,strategyIndex)),mean(frameFeasible(:,strategyIndex)), ...
                        smRatio,alamoutiRatio,qpskRatio,qam16Ratio,metrics.BER, ...
                        metrics.BER_CI95_Lower,metrics.BER_CI95_Upper,bootstrapCI(1),bootstrapCI(2), ...
                        metrics.BLER,mean(frameSE(:,strategyIndex)),mean(frameGoodput(:,strategyIndex)), ...
                        mean(frameNMSE(:,strategyIndex)),mean(frameICI(:,strategyIndex)), ...
                        totalBits,numFrames, ...
                        'VariableNames',{'Channel','SpatialProfile','TxRho','RxRho','MaxDoppler_Hz', ...
                        'PilotSNR_dB','PilotSpacingSymbols','TargetBER','Strategy','MeanPredictedBER', ...
                        'PredictedTargetFeasibleFrameRatio','SpatialMultiplexingFrameRatio', ...
                        'AlamoutiFrameRatio','QPSKFrameRatio','QAM16FrameRatio','BER', ...
                        'BER_WilsonCI95_Lower','BER_WilsonCI95_Upper','BER_FrameBootstrapCI95_Lower', ...
                        'BER_FrameBootstrapCI95_Upper','FER','NetSpectralEfficiency_bpsHz', ...
                        'UncodedGoodput_bpsHz','MeanChannelNMSE','MeanICI_Ratio','PayloadBits','NumFrames'});
                    if isempty(results), results = row; else, results = [results;row]; end %#ok<AGROW>
                    fprintf('  %-22s BER=%8.4g | FER=%5.3f | net SE=%5.2f | feasible=%4.0f%%\n', ...
                        strategies{strategyIndex},metrics.BER,metrics.BLER, ...
                        mean(frameSE(:,strategyIndex)),100*mean(frameFeasible(:,strategyIndex)));
                end
                fprintf('  Adaptive mix: SM %.0f%% / Alamouti %.0f%% | QPSK %.0f%% / 16-QAM %.0f%%\n', ...
                    100*mean(modeIndex==1),100*mean(modeIndex==2), ...
                    100*mean(modulationOrder==4),100*mean(modulationOrder==16));
            end
        end
    end
end

writetable(results,fullfile(outDir,'adaptive_link_results.csv'));
writetable(frameLog,fullfile(outDir,'adaptive_link_frame_log.csv'));
save(fullfile(outDir,'adaptive_link_results.mat'), ...
    'cfg','channelNames','profileNames','profileRhos','dopplerValuesHz', ...
    'pilotSNRValuesdB','pilotSpacing','nSymbols','modulations','targetBER', ...
    'numFrames','numBootstrap','confidenceLevel','strategies','results','frameLog','-v7.3');
plot_pilot_mapped_adaptive_results(results,outDir);
fprintf('\nV15 PILOT-MAPPED ADAPTIVE VALIDATION COMPLETED\n');
fprintf('Channels               : %s\n',strjoin(channelNames,', '));
fprintf('Spatial profiles       : %d\n',numel(profileNames));
fprintf('Doppler points         : %d\n',numel(dopplerValuesHz));
fprintf('Pilot/data SNR points  : %s dB\n',mat2str(pilotSNRValuesdB));
fprintf('Frames/condition       : %d\n',numFrames);
fprintf('Summary rows           : %d (expected %d)\n',height(results),conditionCount*numel(strategies));
fprintf('Frame log rows         : %d (expected %d)\n',height(frameLog),conditionCount*numFrames);
fprintf('Output                 : %s\n',outDir);
end
