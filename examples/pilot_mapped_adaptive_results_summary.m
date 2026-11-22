function summary = pilot_mapped_adaptive_results_summary(results)
%PILOT_MAPPED_ADAPTIVE_RESULTS_SUMMARY Print the main V15 link-adaptation findings.
if nargin < 1 || isempty(results)
    error('pilot_mapped_adaptive_results_summary:MissingResults', ...
        'Pass the table returned by pilot_mapped_adaptive_validation_run.');
end
required = {'Channel','SpatialProfile','MaxDoppler_Hz','PilotSNR_dB','Strategy', ...
    'BER','FER','NetSpectralEfficiency_bpsHz','UncodedGoodput_bpsHz', ...
    'SpatialMultiplexingFrameRatio','AlamoutiFrameRatio','QPSKFrameRatio', ...
    'QAM16FrameRatio','PredictedTargetFeasibleFrameRatio','MeanPredictedBER', ...
    'BER_WilsonCI95_Upper','BER_FrameBootstrapCI95_Upper'};
if ~istable(results) || ~all(ismember(required,results.Properties.VariableNames))
    error('pilot_mapped_adaptive_results_summary:InvalidTable', ...
        'Input table does not contain the required V15 result columns.');
end
fprintf('\n============================================================\n');
fprintf('V15 PILOT-MAPPED ADAPTIVE MIMO-OFDM RESULTS: %d rows\n',height(results));
fprintf('============================================================\n');

channels = unique(string(results.Channel),'stable');
strategies = unique(string(results.Strategy),'stable');
summary = [];

fprintf('\n--- A. Overall results by channel and strategy ---\n');
for c = 1:numel(channels)
    for s = 1:numel(strategies)
        mask = string(results.Channel)==channels(c) & string(results.Strategy)==strategies(s);
        if ~any(mask), continue; end
        meanBER = mean(results.BER(mask));
        meanFER = mean(results.FER(mask));
        meanSE = mean(results.NetSpectralEfficiency_bpsHz(mask));
        meanGoodput = mean(results.UncodedGoodput_bpsHz(mask));
        meanPredictedBER = mean(results.MeanPredictedBER(mask));
        meanFeasible = mean(results.PredictedTargetFeasibleFrameRatio(mask));
        maxWilsonUpper = max(results.BER_WilsonCI95_Upper(mask));
        maxBootstrapUpper = max(results.BER_FrameBootstrapCI95_Upper(mask));
        fprintf(['%-8s | %-22s | BER=%8.4g | FER=%5.3f | SE=%5.2f | goodput=%5.2f ' ...
            '| predBER=%8.4g | predFeas=%4.0f%% | maxU95(W/cluster)=%7.4g/%7.4g\n'], ...
            char(channels(c)),char(strategies(s)),meanBER,meanFER,meanSE,meanGoodput, ...
            meanPredictedBER,100*meanFeasible,maxWilsonUpper,maxBootstrapUpper);
        row = table(channels(c),strategies(s),meanBER,meanFER,meanSE,meanGoodput, ...
            meanPredictedBER,meanFeasible,maxWilsonUpper,maxBootstrapUpper, ...
            'VariableNames',{'Channel','Strategy','MeanBER','MeanFER','MeanNetSE_bpsHz', ...
            'MeanUncodedGoodput_bpsHz','MeanPredictedBER','MeanPredictedFeasibleRatio', ...
            'MaxWilsonUpper95','MaxFrameBootstrapUpper95'});
        if isempty(summary), summary=row; else, summary=[summary;row]; end %#ok<AGROW>
    end
end

fprintf('\n--- B. Adaptive selection mix by channel and pilot/data SNR ---\n');
for c = 1:numel(channels)
    snrs = unique(results.PilotSNR_dB(string(results.Channel)==channels(c)));
    for q = 1:numel(snrs)
        mask = string(results.Channel)==channels(c) & string(results.Strategy)=='adaptive' & ...
            results.PilotSNR_dB==snrs(q);
        if ~any(mask), continue; end
        sm = mean(results.SpatialMultiplexingFrameRatio(mask));
        al = mean(results.AlamoutiFrameRatio(mask));
        qpsk = mean(results.QPSKFrameRatio(mask));
        qam16Values = results.QAM16FrameRatio;
        qam16 = mean(qam16Values(mask));
        feasible = mean(results.PredictedTargetFeasibleFrameRatio(mask));
        ber = mean(results.BER(mask));
        se = mean(results.NetSpectralEfficiency_bpsHz(mask));
        maxWilsonUpper = max(results.BER_WilsonCI95_Upper(mask));
        maxBootstrapUpper = max(results.BER_FrameBootstrapCI95_Upper(mask));
        fprintf(['%-8s | SNR=%2g dB | SM=%4.0f%% Alamouti=%4.0f%% | QPSK=%4.0f%% 16QAM=%4.0f%% ' ...
            '| feasible=%4.0f%% | BER=%8.4g | SE=%5.2f | maxU95(W/cluster)=%7.4g/%7.4g\n'], ...
            char(channels(c)),snrs(q),100*sm,100*al,100*qpsk,100*qam16,100*feasible,ber,se, ...
            maxWilsonUpper,maxBootstrapUpper);
    end
end

fprintf('\n--- C. Adaptive link relative to fixed SM-16QAM ---\n');
for c = 1:numel(channels)
    snrs = unique(results.PilotSNR_dB(string(results.Channel)==channels(c)));
    for q = 1:numel(snrs)
        adaptiveMask = string(results.Channel)==channels(c) & ...
            string(results.Strategy)=='adaptive' & results.PilotSNR_dB==snrs(q);
        fixedMask = string(results.Channel)==channels(c) & ...
            string(results.Strategy)=='fixed-sm-16qam' & results.PilotSNR_dB==snrs(q);
        if ~any(adaptiveMask) || ~any(fixedMask), continue; end
        adaptiveBER = mean(results.BER(adaptiveMask));
        fixedBER = mean(results.BER(fixedMask));
        adaptiveSE = mean(results.NetSpectralEfficiency_bpsHz(adaptiveMask));
        fixedSE = mean(results.NetSpectralEfficiency_bpsHz(fixedMask));
        berGain = 100*(fixedBER-adaptiveBER)/max(fixedBER,eps);
        fprintf('%-8s | SNR=%2g dB | BER change vs fixed SM16=%+7.2f%% | SE adaptive/fixed=%.2f/%.2f\n', ...
            char(channels(c)),snrs(q),berGain,adaptiveSE,fixedSE);
    end
end
fprintf(['\nNote: BER is uncoded. For each channel/SNR group, maxU95 is the largest upper ' ...
    'bound among its profile/Doppler conditions, not a pooled confidence interval.\n']);
fprintf('Predicted feasibility is a link-adaptation estimate, not proof of target-BER compliance.\n');
end
