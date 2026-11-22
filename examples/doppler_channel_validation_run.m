function [results, acfResults] = doppler_channel_validation_run()
%DOPPLER_CHANNEL_VALIDATION_RUN Validate sample-varying Rayleigh/Rician MIMO channels.
% Runs temporal-correlation and within-symbol ICI studies for several Doppler
% frequencies. Output is saved under doppler_validation_results/.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(projectRoot));
cfg = default_config();
cfg.system.subcarrierSpacingHz = 15e3;
cfg.channel.doppler.numOscillators = 64;
cfg.channel.doppler.losDopplerFraction = 1.0;
channelNames = {'rayleigh','rician'};
dopplerValuesHz = [0 30 120 300];
numOFDMSymbols = 192;
samplesPerSymbol = cfg.system.nFFT+cfg.system.cpLength;
sampleRateHz = cfg.system.nFFT*cfg.system.subcarrierSpacingHz;
totalSamples = samplesPerSymbol*numOFDMSymbols;
lagSymbols = [0 1 4 16 64 128];
lagSamples = lagSymbols*samplesPerSymbol;
activeTone = 8;

outDir = fullfile(projectRoot,'doppler_validation_results');
if ~exist(outDir,'dir'), mkdir(outDir); end
results = [];
acfResults = [];

for c = 1:numel(channelNames)
    channelName = channelNames{c};
    for d = 1:numel(dopplerValuesHz)
        fd = dopplerValuesHz(d);
        cfg.channel.doppler.maxDopplerHz = fd;
        seed = cfg.simulation.randomSeed + 91000 + 100*c + d;
        [tapSequence,components,meta] = generate_time_varying_mimo_channel( ...
            cfg,channelName,totalSamples,sampleRateHz,seed);

        [estimatedACF,theoreticalACF] = doppler_temporal_metrics( ...
            components.diffuse,fd,sampleRateHz,lagSamples);
        acfPart = table(repmat(string(channelName),numel(lagSymbols),1), ...
            repmat(fd,numel(lagSymbols),1),lagSymbols(:),real(estimatedACF), ...
            imag(estimatedACF),theoreticalACF(:),abs(estimatedACF-theoreticalACF(:)), ...
            'VariableNames',{'Channel','MaxDoppler_Hz','LagSymbols', ...
            'EstimatedACF_Real','EstimatedACF_Imag','TheoreticalJakesACF','AbsoluteACFError'});
        if isempty(acfResults), acfResults = acfPart; else, acfResults = [acfResults;acfPart]; end %#ok<AGROW>

        % Single active subcarrier probe: nonzero off-tone energy is ICI.
        Xprobe = complex(zeros(cfg.system.nFFT,cfg.system.nTx,numOFDMSymbols));
        Xprobe(activeTone,1,:) = 1;
        txWaveform = ofdm_modulate(Xprobe,cfg.system.cpLength);
        rxWaveform = apply_time_varying_mimo_ofdm_channel( ...
            txWaveform,tapSequence,0);
        Yprobe = ofdm_demodulate(rxWaveform,cfg.system.nFFT,cfg.system.cpLength);
        totalPower = sum(abs(Yprobe(:)).^2);
        desiredPower = sum(sum(abs(reshape(Yprobe(activeTone,:,:), ...
            cfg.system.nRx,numOFDMSymbols)).^2));
        iciRatio = max(totalPower-desiredPower,0)/max(totalPower,eps);

        diffusePower = mean(abs(components.diffuse(:)).^2);
        losPower = mean(abs(components.los(:)).^2);
        acfError = mean(abs(estimatedACF-theoreticalACF(:)));
        row = table(string(channelName),fd,fd/cfg.system.subcarrierSpacingHz, ...
            meta.coherenceTimeApprox_s,sampleRateHz,numOFDMSymbols, ...
            acfError,iciRatio,diffusePower,losPower, ...
            'VariableNames',{'Channel','MaxDoppler_Hz','NormalizedDoppler', ...
            'ApproxCoherenceTime_s','SampleRate_Hz','NumOFDMSymbols', ...
            'MeanAbsoluteACFError','ICI_LeakageRatio','MeanDiffusePower', ...
            'MeanLOSPower'});
        if isempty(results), results = row; else, results = [results;row]; end %#ok<AGROW>

        fprintf(['%s | fD=%3g Hz | ACF MAE=%.4g | ICI=%.4g | ' ...
            'Tc~%.4g s\n'],upper(channelName),fd,acfError,iciRatio,meta.coherenceTimeApprox_s);
    end
end

writetable(results,fullfile(outDir,'doppler_channel_validation.csv'));
writetable(acfResults,fullfile(outDir,'doppler_temporal_acf.csv'));
save(fullfile(outDir,'doppler_channel_validation.mat'), ...
    'cfg','channelNames','dopplerValuesHz','results','acfResults');
plot_doppler_validation(results,acfResults,outDir);
fprintf('\nDOPPLER CHANNEL VALIDATION COMPLETED\n');
fprintf('Channels             : %s\n',strjoin(channelNames,', '));
fprintf('Doppler points       : %d\n',numel(dopplerValuesHz));
fprintf('OFDM symbols/run     : %d\n',numOFDMSymbols);
fprintf('Result rows          : %d\n',height(results));
fprintf('Output               : %s\n',outDir);
end
