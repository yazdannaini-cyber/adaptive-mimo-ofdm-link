function results = pilot_mapped_adaptive_validation_run_v15_1()
%PILOT_MAPPED_ADAPTIVE_VALIDATION_RUN_V15_1 V15.1 burst-adaptive link study.
% Decisions are refreshed after every orthogonal pilot occasion. MCS metrics
% include temporal-LMMSE posterior channel uncertainty and causal pilot-ICI
% innovation estimates. Output includes summary, frame and block-decision logs.

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
numFrames = 8;       % Screening run; raise to 100+ for confidence-grade studies.
numBootstrap = 2000;
confidenceLevel = 0.95;
strategies = {'adaptive','fixed-sm-16qam','fixed-alamouti-qpsk'};

outDir = fullfile(projectRoot,'pilot_mapped_results_v15_1');
if ~exist(outDir,'dir'), mkdir(outDir); end
results = [];
frameLog = [];
blockLog = [];
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
                fprintf('\n[%d/%d] %s | %s | fD=%g Hz | SNR=%g dB\n', ...
                    conditionNumber,conditionCount,upper(channelName),profileName,fd,snrDb);

                frameErrors = zeros(numFrames,numel(strategies));
                frameBits = zeros(numFrames,numel(strategies));
                frameFER = false(numFrames,numel(strategies));
                frameSE = zeros(numFrames,numel(strategies));
                frameGoodput = zeros(numFrames,numel(strategies));
                frameNMSE = zeros(numFrames,numel(strategies));
                frameICI = zeros(numFrames,numel(strategies));
                framePredictedBER = zeros(numFrames,numel(strategies));
                frameFeasibleRatio = zeros(numFrames,numel(strategies));
                frameSMRatio = zeros(numFrames,numel(strategies));
                frameAlamoutiRatio = zeros(numFrames,numel(strategies));
                frameQPSKRatio = zeros(numFrames,numel(strategies));
                frameQAM16Ratio = zeros(numFrames,numel(strategies));
                framePosterior = zeros(numFrames,numel(strategies));
                frameEffectiveNoise = zeros(numFrames,numel(strategies));
                frameEstimatedICI = zeros(numFrames,numel(strategies));

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

                    frame = simulate_pilot_mapped_adaptive_frame_v15_1(cfg,channelName,params);
                    candidates = {frame.adaptive,frame.fixedSM16QAM,frame.fixedAlamoutiQPSK};
                    for strategyIndex = 1:numel(strategies)
                        candidate = candidates{strategyIndex};
                        frameErrors(frameIndex,strategyIndex) = candidate.BitErrors;
                        frameBits(frameIndex,strategyIndex) = candidate.TotalBits;
                        frameFER(frameIndex,strategyIndex) = candidate.FrameError;
                        frameSE(frameIndex,strategyIndex) = candidate.NetSpectralEfficiency_bpsHz;
                        frameGoodput(frameIndex,strategyIndex) = candidate.UncodedGoodput_bpsHz;
                        frameNMSE(frameIndex,strategyIndex) = candidate.NMSE;
                        frameICI(frameIndex,strategyIndex) = candidate.ICI_Ratio;
                        framePredictedBER(frameIndex,strategyIndex) = candidate.PredictedBER;
                        frameFeasibleRatio(frameIndex,strategyIndex) = candidate.TargetFeasibleBlockRatio;
                        frameSMRatio(frameIndex,strategyIndex) = candidate.SpatialMultiplexingFrameRatio;
                        frameAlamoutiRatio(frameIndex,strategyIndex) = candidate.AlamoutiFrameRatio;
                        frameQPSKRatio(frameIndex,strategyIndex) = candidate.QPSKFrameRatio;
                        frameQAM16Ratio(frameIndex,strategyIndex) = candidate.QAM16FrameRatio;
                        framePosterior(frameIndex,strategyIndex) = candidate.MeanPosteriorVariance;
                        frameEffectiveNoise(frameIndex,strategyIndex) = candidate.MeanEffectiveNoiseVariance;
                        frameEstimatedICI(frameIndex,strategyIndex) = candidate.MeanEstimatedICINoiseVariance;
                    end

                    logRow = table(string(channelName),string(profileName),fd,snrDb, ...
                        frameIndex,string(frame.adaptive.Mode),logical(frame.adaptive.CausalReceiverEstimates), ...
                        frame.adaptive.ModulationOrder, ...
                        frame.adaptive.PredictedBER,frame.adaptive.TargetFeasibleBlockRatio, ...
                        frame.adaptive.SpatialMultiplexingFrameRatio,frame.adaptive.AlamoutiFrameRatio, ...
                        frame.adaptive.QPSKFrameRatio,frame.adaptive.QAM16FrameRatio, ...
                        frame.adaptive.BER,frame.adaptive.NetSpectralEfficiency_bpsHz, ...
                        frame.adaptive.NMSE,frame.adaptive.MeanPosteriorVariance, ...
                        frame.adaptive.MeanEffectiveNoiseVariance,frame.adaptive.MeanEstimatedICINoiseVariance, ...
                        frame.adaptive.ICI_Ratio, ...
                        'VariableNames',{'Channel','SpatialProfile','MaxDoppler_Hz','PilotSNR_dB', ...
                        'FrameIndex','SelectedMode','ReceiverUsesCausalCSI','SingleModulationOrder','PredictedBER', ...
                        'PredictedFeasibleDataRatio','SMDataRatio','AlamoutiDataRatio', ...
                        'QPSKDataRatio','QAM16DataRatio','AdaptiveBER','AdaptiveNetSE_bpsHz', ...
                        'AdaptiveNMSE','MeanPosteriorVariance','MeanEffectiveNoiseVariance', ...
                        'MeanEstimatedICINoiseVariance','AdaptiveICI_Ratio'});
                    if isempty(frameLog), frameLog = logRow; else, frameLog = [frameLog;logRow]; end %#ok<AGROW>

                    for b = 1:numel(frame.blocks)
                        block = frame.blocks(b);
                        blockChoice = block.choice;
                        bRow = table(string(channelName),string(profileName),fd,snrDb,frameIndex,b, ...
                            block.dataIndices(1),block.dataIndices(end),string(blockChoice.Mode), ...
                            blockChoice.ModulationOrder,blockChoice.PredictedBER,blockChoice.TargetFeasible, ...
                            block.prediction.predictedBER(1,1),block.prediction.predictedBER(1,2), ...
                            block.prediction.predictedBER(2,1),block.prediction.predictedBER(2,2), ...
                            block.prediction.meanPosteriorVariance,block.iciNoiseVariance, ...
                            block.iciDetail.analyticalICIFloor, ...
                            'VariableNames',{'Channel','SpatialProfile','MaxDoppler_Hz','PilotSNR_dB', ...
                            'FrameIndex','BurstIndex','FirstDataSymbol','LastDataSymbol','SelectedMode', ...
                            'SelectedModulationOrder','PredictedBER','TargetFeasible', ...
                            'PredBER_SM_QPSK','PredBER_SM_16QAM','PredBER_Alamouti_QPSK', ...
                            'PredBER_Alamouti_16QAM','MeanPosteriorVariance','EstimatedICINoiseVariance', ...
                            'AnalyticalICIFloor'});
                        if isempty(blockLog), blockLog=bRow; else, blockLog=[blockLog;bRow]; end %#ok<AGROW>
                    end
                end

                for strategyIndex = 1:numel(strategies)
                    totalErrors = sum(frameErrors(:,strategyIndex));
                    totalBits = sum(frameBits(:,strategyIndex));
                    totalFrameErrors = nnz(frameFER(:,strategyIndex));
                    metrics = bit_error_metrics(totalErrors,totalBits,totalFrameErrors,numFrames);
                    bootstrapCI = cluster_bootstrap_ber_ci(frameErrors(:,strategyIndex), ...
                        frameBits(:,strategyIndex),numBootstrap,confidenceLevel, ...
                        cfg.statistics.bootstrapSeed+conditionNumber*100+strategyIndex);
                    row = table(string(channelName),string(profileName),rhoTx,rhoRx,fd,snrDb, ...
                        pilotSpacing,targetBER,string(strategies{strategyIndex}), ...
                        mean(framePredictedBER(:,strategyIndex)),mean(frameFeasibleRatio(:,strategyIndex)), ...
                        mean(frameSMRatio(:,strategyIndex)),mean(frameAlamoutiRatio(:,strategyIndex)), ...
                        mean(frameQPSKRatio(:,strategyIndex)),mean(frameQAM16Ratio(:,strategyIndex)), ...
                        metrics.BER,metrics.BER_CI95_Lower,metrics.BER_CI95_Upper,bootstrapCI(1),bootstrapCI(2), ...
                        metrics.BLER,mean(frameSE(:,strategyIndex)),mean(frameGoodput(:,strategyIndex)), ...
                        mean(frameNMSE(:,strategyIndex)),mean(framePosterior(:,strategyIndex)), ...
                        mean(frameEffectiveNoise(:,strategyIndex)),mean(frameEstimatedICI(:,strategyIndex)), ...
                        mean(frameICI(:,strategyIndex)),totalBits,numFrames, ...
                        'VariableNames',{'Channel','SpatialProfile','TxRho','RxRho','MaxDoppler_Hz', ...
                        'PilotSNR_dB','PilotSpacingSymbols','TargetBER','Strategy','MeanPredictedBER', ...
                        'PredictedTargetFeasibleDataRatio','SpatialMultiplexingDataRatio', ...
                        'AlamoutiDataRatio','QPSKDataRatio','QAM16DataRatio','BER', ...
                        'BER_WilsonCI95_Lower','BER_WilsonCI95_Upper','BER_FrameBootstrapCI95_Lower', ...
                        'BER_FrameBootstrapCI95_Upper','FER','NetSpectralEfficiency_bpsHz', ...
                        'UncodedGoodput_bpsHz','MeanChannelNMSE','MeanPosteriorVariance', ...
                        'MeanEffectiveNoiseVariance','MeanEstimatedICINoiseVariance','MeanICI_Ratio', ...
                        'PayloadBits','NumFrames'});
                    if isempty(results), results=row; else, results=[results;row]; end %#ok<AGROW>
                    fprintf('  %-22s BER=%8.4g | FER=%5.3f | netSE=%5.2f | predFeas=%4.0f%%\n', ...
                        strategies{strategyIndex},metrics.BER,metrics.BLER, ...
                        mean(frameSE(:,strategyIndex)),100*mean(frameFeasibleRatio(:,strategyIndex)));
                end
                fprintf('  Burst-adaptive mix: SM %.0f%% / Alamouti %.0f%% | QPSK %.0f%% / 16-QAM %.0f%%\n', ...
                    100*mean(frameSMRatio(:,1)),100*mean(frameAlamoutiRatio(:,1)), ...
                    100*mean(frameQPSKRatio(:,1)),100*mean(frameQAM16Ratio(:,1)));
            end
        end
    end
end

writetable(results,fullfile(outDir,'adaptive_link_results_v15_1.csv'));
writetable(frameLog,fullfile(outDir,'adaptive_link_frame_log_v15_1.csv'));
writetable(blockLog,fullfile(outDir,'adaptive_link_block_decisions_v15_1.csv'));
save(fullfile(outDir,'adaptive_link_results_v15_1.mat'), ...
    'cfg','channelNames','profileNames','profileRhos','dopplerValuesHz', ...
    'pilotSNRValuesdB','pilotSpacing','nSymbols','modulations','targetBER', ...
    'numFrames','numBootstrap','confidenceLevel','strategies','results','frameLog','blockLog','-v7.3');
plot_pilot_mapped_adaptive_results_v15_1(results,outDir);
fprintf('\nV15.1 PILOT-MAPPED BURST-ADAPTIVE VALIDATION COMPLETED\n');
fprintf('Channels               : %s\n',strjoin(channelNames,', '));
fprintf('Spatial profiles       : %d\n',numel(profileNames));
fprintf('Doppler points         : %d\n',numel(dopplerValuesHz));
fprintf('Pilot/data SNR points  : %s dB\n',mat2str(pilotSNRValuesdB));
fprintf('Frames/condition       : %d\n',numFrames);
fprintf('Summary rows           : %d (expected %d)\n',height(results),conditionCount*numel(strategies));
fprintf('Frame log rows         : %d (expected %d)\n',height(frameLog),conditionCount*numFrames);
pilotScheduleForCount = mimo_ofdm_pilot_schedule(nSymbols,cfg.system.nTx,pilotSpacing);
fprintf('Block decision rows    : %d (expected %d)\n',height(blockLog),conditionCount*numFrames*numel(pilotScheduleForCount.pilotOccasionStarts));
fprintf('Output                 : %s\n',outDir);
end
