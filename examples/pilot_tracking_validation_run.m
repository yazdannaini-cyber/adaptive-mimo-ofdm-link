function results = pilot_tracking_validation_run()
%PILOT_TRACKING_VALIDATION_RUN Benchmark temporal pilot tracking under Doppler.
% The experiment compares causal hold-last, noncausal linear interpolation,
% and a causal Jakes-covariance LMMSE predictor. Channel snapshots are sampled
% once per OFDM-symbol period; intra-symbol ICI remains a separate V11 study.
%
% Each pilot epoch provides orthogonal per-link LS observations on the full
% frequency grid. This isolates temporal tracking from pilot waveform design.
% The known Rician LOS trajectory is supplied as a model prior to all methods.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
cfg = default_config();
channelNames = {'rayleigh','rician'};
dopplerValuesHz = [0 30 120 300];
profileNames = {'iid','joint-moderate','joint-high'};
profileRhos = [0 0; 0.7 0.7; 0.9 0.9];
pilotSpacingValues = [2 4 8 16];
pilotSNRValuesdB = [0 10 20];
methods = {'hold-last','linear','jakes-lmmse'};
nSymbols = 64;
numFrames = 40;
numBootstrap = 1000;
confidenceLevel = 0.95;
N = cfg.system.nFFT;
Nr = cfg.system.nRx;
Nt = cfg.system.nTx;
Tsym = (N+cfg.system.cpLength)/(N*cfg.system.subcarrierSpacingHz);
symbolRateHz = 1/Tsym;

outDir = fullfile(projectRoot,'pilot_tracking_results');
if ~exist(outDir,'dir'), mkdir(outDir); end
results = [];
frameRecords = struct();
rowCount = 0;
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>

for c = 1:numel(channelNames)
    channelName = channelNames{c};
    if strcmp(channelName,'rayleigh')
        diffuseVariance = 1;
        meanChannelPower = 1;
    else
        K = 10^(cfg.channel.ricianK_dB/10);
        diffuseVariance = 1/(K+1);
        % The configured LOS component occupies the first PDP tap only.
        meanChannelPower = diffuseVariance + cfg.channel.pdp(1)*K/(K+1);
    end

    for profile = 1:numel(profileNames)
        cfgProfile = cfg;
        cfgProfile.channel.spatialCorrelation.enabled = true;
        cfgProfile.channel.spatialCorrelation.txRho = profileRhos(profile,1);
        cfgProfile.channel.spatialCorrelation.rxRho = profileRhos(profile,2);

        for d = 1:numel(dopplerValuesHz)
            fd = dopplerValuesHz(d);
            cfgFrame = cfgProfile;
            cfgFrame.channel.doppler.maxDopplerHz = fd;
            nSpacing = numel(pilotSpacingValues);
            nSNR = numel(pilotSNRValuesdB);
            nMethods = numel(methods);
            frameNmse = nan(numFrames,nSpacing,nSNR,nMethods);
            frameMse = nan(numFrames,nSpacing,nSNR,nMethods);
            framePosteriorMSE = nan(numFrames,nSpacing,nSNR,nMethods);
            pilotIndexCache = cell(nSpacing,1);
            dataIndexCache = cell(nSpacing,1);
            setupCache = cell(nSpacing,nSNR);

            for ip = 1:nSpacing
                spacing = pilotSpacingValues(ip);
                pilotIndices = 1:spacing:nSymbols;
                dataIndices = setdiff(1:nSymbols,pilotIndices);
                pilotIndexCache{ip} = pilotIndices;
                dataIndexCache{ip} = dataIndices;
                for is = 1:nSNR
                    pilotSNRdB = pilotSNRValuesdB(is);
                    noiseVariance = meanChannelPower*10^(-pilotSNRdB/10);
                    setupCache{ip,is} = temporal_channel_tracker_setup( ...
                        pilotIndices,nSymbols,noiseVariance,fd,Tsym,diffuseVariance);
                end
            end

            for frame = 1:numFrames
                channelSeed = cfg.simulation.randomSeed + 2e6*c + ...
                    2e5*profile + 2e4*d + frame;
                [tapSequence,components] = generate_time_varying_mimo_channel( ...
                    cfgFrame,channelName,nSymbols,symbolRateHz,channelSeed);
                Htemp = fft(tapSequence,N,2);
                Htrue = permute(Htemp,[2 3 4 1]);
                HmeanTemp = fft(components.los,N,2);
                Hmean = permute(HmeanTemp,[2 3 4 1]);

                for ip = 1:nSpacing
                    pilotIndices = pilotIndexCache{ip};
                    dataIndices = dataIndexCache{ip};
                    pilotTruth = Htrue(:,:,:,pilotIndices);
                    for is = 1:nSNR
                        pilotSNRdB = pilotSNRValuesdB(is);
                        noiseVariance = meanChannelPower*10^(-pilotSNRdB/10);
                        noiseSeed = cfg.simulation.randomSeed + 5e6 + ...
                            1e6*c + 1e5*profile + 1e4*d + ...
                            1e3*ip + 100*is + frame;
                        rng(noiseSeed,'twister');
                        pilotNoise = sqrt(noiseVariance/2)*( ...
                            randn(size(pilotTruth))+1i*randn(size(pilotTruth)));
                        pilotObservations = pilotTruth+pilotNoise;

                        [Hhold,~] = temporal_channel_tracking_apply( ...
                            pilotObservations,pilotIndices,nSymbols,'hold-last',Hmean,[]);
                        [Hlinear,~] = temporal_channel_tracking_apply( ...
                            pilotObservations,pilotIndices,nSymbols,'linear',Hmean,[]);
                        setup = setupCache{ip,is};
                        [Hjakes,posteriorVariance] = temporal_channel_tracking_apply( ...
                            pilotObservations,pilotIndices,nSymbols,'jakes-lmmse',Hmean,setup);
                        estimates = {Hhold,Hlinear,Hjakes};
                        for im = 1:nMethods
                            estimate = estimates{im}(:,:,:,dataIndices);
                            truth = Htrue(:,:,:,dataIndices);
                            errorArray = estimate-truth;
                            frameMse(frame,ip,is,im) = mean(abs(errorArray(:)).^2);
                            frameNmse(frame,ip,is,im) = sum(abs(errorArray(:)).^2)/ ...
                                max(sum(abs(truth(:)).^2),eps);
                            if im == 3
                                framePosteriorMSE(frame,ip,is,im) = ...
                                    mean(posteriorVariance(dataIndices));
                            end
                        end
                    end
                end
            end

            % Reduce independent frame clusters and preserve paired method differences.
            for ip = 1:nSpacing
                pilotIndices = pilotIndexCache{ip};
                dataIndices = dataIndexCache{ip};
                for is = 1:nSNR
                    pilotSNRdB = pilotSNRValuesdB(is);
                    referenceNmse = frameNmse(:,ip,is,1);
                    for im = 1:nMethods
                        valuesNmse = frameNmse(:,ip,is,im);
                        valuesMse = frameMse(:,ip,is,im);
                        rowSeed = cfg.statistics.bootstrapSeed + 1e6*c + ...
                            1e5*profile + 1e4*d + 100*ip + 10*is + im;
                        nmseCI = frame_bootstrap_mean_ci(valuesNmse, ...
                            numBootstrap,confidenceLevel,rowSeed);
                        if im == 1
                            improvementPercent = 0;
                            pairedCI = [NaN NaN];
                        else
                            difference = referenceNmse-valuesNmse;
                            pairedCI = frame_bootstrap_mean_ci(difference, ...
                                numBootstrap,confidenceLevel,rowSeed+50000);
                            improvementPercent = 100*(mean(referenceNmse)-mean(valuesNmse))/ ...
                                max(mean(referenceNmse),eps);
                        end
                        if im == 3
                            posteriorMSE = mean(framePosteriorMSE(:,ip,is,im));
                            calibrationRatio = mean(valuesMse)/max(posteriorMSE,eps);
                        else
                            posteriorMSE = NaN;
                            calibrationRatio = NaN;
                        end
                        row = table(string(channelName),string(profileNames{profile}), ...
                            profileRhos(profile,1),profileRhos(profile,2),fd, ...
                            pilotSpacingValues(ip),pilotSNRdB,string(methods{im}), ...
                            numFrames,nSymbols,numel(pilotIndices),numel(dataIndices), ...
                            mean(valuesNmse),10*log10(max(mean(valuesNmse),realmin)), ...
                            nmseCI(1),nmseCI(2),mean(valuesMse),posteriorMSE, ...
                            calibrationRatio,improvementPercent,pairedCI(1),pairedCI(2), ...
                            'VariableNames',{'Channel','SpatialProfile','TxRho','RxRho', ...
                            'MaxDoppler_Hz','PilotSpacingSymbols','PilotSNR_dB','Estimator', ...
                            'NumFrames','NumSymbols','NumPilotSymbols','NumDataSymbols', ...
                            'MeanNMSE','MeanNMSE_dB','NMSE_CI95_Lower','NMSE_CI95_Upper', ...
                            'MeanMSE_PerCoefficient','MeanPosteriorMSE_PerCoefficient', ...
                            'PosteriorMSE_CalibrationRatio','NMSEImprovementVsHold_Percent', ...
                            'PairedNMSEDifferenceVsHold_CI95_Lower', ...
                            'PairedNMSEDifferenceVsHold_CI95_Upper'});
                        if isempty(results), results = row; else, results = [results;row]; end %#ok<AGROW>
                        rowCount = rowCount+1;
                    end
                end
            end

            key = sprintf('c%d_p%d_d%d',c,profile,d);
            frameRecords.(key).nmse = frameNmse;
            frameRecords.(key).mse = frameMse;
            frameRecords.(key).posteriorMSE = framePosteriorMSE;
        end
    end
end

writetable(results,fullfile(outDir,'pilot_tracking_validation.csv'));
save(fullfile(outDir,'pilot_tracking_validation.mat'), ...
    'cfg','channelNames','dopplerValuesHz','profileNames','profileRhos', ...
    'pilotSpacingValues','pilotSNRValuesdB','methods','nSymbols','numFrames', ...
    'numBootstrap','confidenceLevel','Tsym','results','frameRecords','-v7.3');
plot_pilot_tracking_validation(results,outDir);
fprintf('\nPILOT TRACKING VALIDATION COMPLETED\n');
fprintf('Channels              : %s\n',strjoin(channelNames,', '));
fprintf('Spatial profiles     : %d\n',numel(profileNames));
fprintf('Doppler points        : %d\n',numel(dopplerValuesHz));
fprintf('Pilot spacings        : %s symbols\n',mat2str(pilotSpacingValues));
fprintf('Pilot SNR points      : %s dB\n',mat2str(pilotSNRValuesdB));
fprintf('Frames/condition      : %d\n',numFrames);
fprintf('Estimators             : %s\n',strjoin(methods,', '));
fprintf('Result rows            : %d\n',rowCount);
fprintf('Output                 : %s\n',outDir);
end
