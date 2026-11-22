function results = pilot_mapped_ofdm_validation_run()
%PILOT_MAPPED_OFDM_VALIDATION_RUN End-to-end pilot-mapped OFDM study under Doppler.
% Unlike V12's isolated LS-observation study, this experiment maps explicit
% orthogonal pilot OFDM symbols and QAM payloads, applies a sample-varying MIMO
% FIR channel, demodulates the received waveform, tracks the channel from the
% extracted pilots, and evaluates data-symbol NMSE/BER. ICI is measured from
% the difference between the actual waveform output and the diagonal-channel
% reference; it is not yet cancelled by this V13 receiver.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
cfg = default_config();
cfg.channel.doppler.numOscillators = 32;

channelNames = {'rayleigh','rician'};
profileNames = {'iid','joint-high'};
profileRhos = [0 0; 0.9 0.9];
dopplerValuesHz = [0 120 300];
pilotSpacingValues = [4 8 16];
pilotSNRdB = 15;
modulationOrder = 16;
nSymbols = 24;
numFrames = 8;
numBootstrap = 500;
confidenceLevel = 0.95;
methods = {'hold-last','linear','jakes-lmmse','perfect-csi','perfect-csi-ici'};

outDir = fullfile(projectRoot,'pilot_mapped_results');
if ~exist(outDir,'dir'), mkdir(outDir); end
results = [];
rowCount = 0;
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>

conditionCount = numel(channelNames)*numel(profileNames)* ...
    numel(dopplerValuesHz)*numel(pilotSpacingValues);
conditionNumber = 0;
for c = 1:numel(channelNames)
    channelName = channelNames{c};
    for p = 1:numel(profileNames)
        for d = 1:numel(dopplerValuesHz)
            for ip = 1:numel(pilotSpacingValues)
                conditionNumber = conditionNumber+1;
                fd = dopplerValuesHz(d);
                spacing = pilotSpacingValues(ip);
                rhoTx = profileRhos(p,1);
                rhoRx = profileRhos(p,2);
                fprintf('\n[%d/%d] %s | profile=%s | fD=%g Hz | pilot spacing=%d | pilot SNR=%g dB\n', ...
                    conditionNumber,conditionCount,upper(channelName),profileNames{p}, ...
                    fd,spacing,pilotSNRdB);

                frameNMSE = nan(numFrames,numel(methods));
                frameMSE = nan(numFrames,numel(methods));
                framePosteriorMSE = nan(numFrames,numel(methods));
                frameBER = nan(numFrames,numel(methods));
                frameErrors = zeros(numFrames,numel(methods));
                totalBitsPerFrame = zeros(numFrames,1);
                frameFailure = false(numFrames,numel(methods));
                frameICI = nan(numFrames,1);
                framePilotContamination = nan(numFrames,1);
                framePilotOverhead = nan(numFrames,1);
                frameSpectralEfficiency = nan(numFrames,1);

                for frameIndex = 1:numFrames
                    seedBase = cfg.simulation.randomSeed + 1e7*c + ...
                        1e6*p + 1e5*d + 1e4*ip + 100*frameIndex;
                    params.nSymbols = nSymbols;
                    params.pilotSpacing = spacing;
                    params.modulationOrder = modulationOrder;
                    params.pilotSNRdB = pilotSNRdB;
                    params.maxDopplerHz = fd;
                    params.txRho = rhoTx;
                    params.rxRho = rhoRx;
                    params.dataSeed = seedBase+1;
                    params.pilotSeed = seedBase+2;
                    params.channelSeed = seedBase+3;
                    params.noiseSeed = seedBase+4;
                    frame = simulate_pilot_mapped_ofdm_frame(cfg,channelName,params);

                    frameICI(frameIndex) = frame.ICI_Ratio;
                    framePilotContamination(frameIndex) = frame.PilotContaminationVariance;
                    framePilotOverhead(frameIndex) = frame.PilotOverhead;
                    frameSpectralEfficiency(frameIndex) = frame.SpectralEfficiency;
                    totalBitsPerFrame(frameIndex) = frame.NumBits;
                    for m = 1:numel(methods)
                        frameNMSE(frameIndex,m) = frame.methods(m).NMSE;
                        frameMSE(frameIndex,m) = frame.methods(m).MSEPerCoefficient;
                        framePosteriorMSE(frameIndex,m) = frame.methods(m).PosteriorMSEPerCoefficient;
                        frameBER(frameIndex,m) = frame.methods(m).BER;
                        frameErrors(frameIndex,m) = frame.methods(m).BitErrors;
                        frameFailure(frameIndex,m) = frame.methods(m).FrameError;
                    end
                end

                for m = 1:numel(methods)
                    totalErrors = sum(frameErrors(:,m));
                    totalBits = sum(totalBitsPerFrame);
                    frameErrorsCount = nnz(frameFailure(:,m));
                    metrics = bit_error_metrics(totalErrors,totalBits,frameErrorsCount,numFrames);
                    nmseCI = frame_bootstrap_mean_ci(frameNMSE(:,m),numBootstrap, ...
                        confidenceLevel,cfg.statistics.bootstrapSeed+conditionNumber*100+m);
                    berClusterCI = frame_bootstrap_mean_ci(frameBER(:,m),numBootstrap, ...
                        confidenceLevel,cfg.statistics.bootstrapSeed+conditionNumber*100+5000+m);
                    meanPosteriorMSE = mean(framePosteriorMSE(:,m),'omitnan');
                    if ~isfinite(meanPosteriorMSE), meanPosteriorMSE = NaN; end
                    calibrationRatio = mean(frameMSE(:,m),'omitnan')/max(meanPosteriorMSE,eps);
                    if ~strcmp(methods{m},'jakes-lmmse')
                        calibrationRatio = NaN;
                        meanPosteriorMSE = NaN;
                    end
                    spectralEfficiency = mean(frameSpectralEfficiency);
                    row = table(string(channelName),string(profileNames{p}),rhoTx,rhoRx, ...
                        fd,spacing,pilotSNRdB,modulationOrder,string(methods{m}), ...
                        numFrames,nSymbols,mean(framePilotOverhead), ...
                        mean(frameICI),mean(framePilotContamination), ...
                        mean(frameNMSE(:,m)),nmseCI(1),nmseCI(2), ...
                        mean(frameMSE(:,m)),meanPosteriorMSE,calibrationRatio, ...
                        metrics.BER,metrics.BER_CI95_Lower,metrics.BER_CI95_Upper, ...
                        berClusterCI(1),berClusterCI(2),metrics.BLER, ...
                        spectralEfficiency,spectralEfficiency*(1-metrics.BER), ...
                        'VariableNames',{'Channel','SpatialProfile','TxRho','RxRho', ...
                        'MaxDoppler_Hz','PilotSpacingSymbols','PilotSNR_dB','ModulationOrder', ...
                        'Estimator','NumFrames','NumOFDMSymbols','PilotOverheadRatio', ...
                        'ICI_Ratio','PilotObservationContaminationVariance', ...
                        'MeanNMSE','NMSE_CI95_Lower','NMSE_CI95_Upper', ...
                        'MeanMSE_PerCoefficient','MeanPosteriorMSE_PerCoefficient', ...
                        'PosteriorMSE_CalibrationRatio','BER','BER_WilsonCI95_Lower', ...
                        'BER_WilsonCI95_Upper','BER_FrameBootstrapCI95_Lower', ...
                        'BER_FrameBootstrapCI95_Upper','FER','NetSpectralEfficiency_bpsHz', ...
                        'UncodedGoodput_bpsHz'});
                    if isempty(results), results = row; else, results = [results;row]; end %#ok<AGROW>
                    rowCount = rowCount+1;
                    fprintf('  %-12s NMSE=%8.4g | BER=%8.4g | FER=%5.3f | net SE=%5.2f | ICI=%7.3g\n', ...
                        methods{m},mean(frameNMSE(:,m)),metrics.BER,metrics.BLER, ...
                        spectralEfficiency,mean(frameICI));
                end
            end
        end
    end
end

writetable(results,fullfile(outDir,'pilot_mapped_ofdm_validation.csv'));
save(fullfile(outDir,'pilot_mapped_ofdm_validation.mat'), ...
    'cfg','channelNames','profileNames','profileRhos','dopplerValuesHz', ...
    'pilotSpacingValues','pilotSNRdB','modulationOrder','nSymbols','numFrames', ...
    'numBootstrap','confidenceLevel','methods','results','-v7.3');
plot_pilot_mapped_ofdm_results(results,outDir);
fprintf('\nPILOT-MAPPED OFDM VALIDATION COMPLETED\n');
fprintf('Channels               : %s\n',strjoin(channelNames,', '));
fprintf('Spatial profiles       : %d\n',numel(profileNames));
fprintf('Doppler points         : %d\n',numel(dopplerValuesHz));
fprintf('Pilot spacings         : %s OFDM symbols\n',mat2str(pilotSpacingValues));
fprintf('Pilot SNR              : %g dB\n',pilotSNRdB);
fprintf('Frames/condition       : %d\n',numFrames);
fprintf('Estimator/mode rows    : %d (expected %d)\n',rowCount, ...
    conditionCount*numel(methods));
fprintf('Output                 : %s\n',outDir);
end
