function summary = pilot_mapped_adaptive_results_summary_v15_1(results)
%PILOT_MAPPED_ADAPTIVE_RESULTS_SUMMARY_V15_1 Summarize burst-adaptive results.
if nargin<1 || isempty(results) || ~istable(results)
    error('pilot_mapped_adaptive_results_summary_v15_1:MissingResults', ...
        'Pass the table returned by pilot_mapped_adaptive_validation_run_v15_1.');
end
required={'Channel','SpatialProfile','MaxDoppler_Hz','PilotSNR_dB','Strategy', ...
    'BER','FER','NetSpectralEfficiency_bpsHz','UncodedGoodput_bpsHz', ...
    'SpatialMultiplexingDataRatio','AlamoutiDataRatio','QPSKDataRatio','QAM16DataRatio', ...
    'PredictedTargetFeasibleDataRatio','MeanPredictedBER','BER_WilsonCI95_Upper', ...
    'BER_FrameBootstrapCI95_Upper','MeanPosteriorVariance','MeanEstimatedICINoiseVariance'};
if ~all(ismember(required,results.Properties.VariableNames))
    error('pilot_mapped_adaptive_results_summary_v15_1:InvalidTable', ...
        'Input table does not contain the required V15.1 result columns.');
end
fprintf('\n============================================================\n');
fprintf('V15.1 BURST-ADAPTIVE PILOT-MAPPED RESULTS: %d rows\n',height(results));
fprintf('============================================================\n');
channels=unique(string(results.Channel),'stable');
strategies=unique(string(results.Strategy),'stable');
summary=[];

fprintf('\n--- A. Overall results by channel and strategy ---\n');
for c=1:numel(channels)
    for s=1:numel(strategies)
        mask=string(results.Channel)==channels(c) & string(results.Strategy)==strategies(s);
        if ~any(mask), continue; end
        meanBER=mean(results.BER(mask));
        meanFER=mean(results.FER(mask));
        meanSE=mean(results.NetSpectralEfficiency_bpsHz(mask));
        meanGoodput=mean(results.UncodedGoodput_bpsHz(mask));
        meanPred=mean(results.MeanPredictedBER(mask));
        feas=mean(results.PredictedTargetFeasibleDataRatio(mask));
        maxW=max(results.BER_WilsonCI95_Upper(mask));
        maxB=max(results.BER_FrameBootstrapCI95_Upper(mask));
        meanPost=mean(results.MeanPosteriorVariance(mask));
        meanICI=mean(results.MeanEstimatedICINoiseVariance(mask));
        fprintf(['%-8s | %-22s | BER=%8.4g | FER=%5.3f | SE=%5.2f | goodput=%5.2f ' ...
            '| predBER=%8.4g | feasible=%4.0f%% | postVar=%7.4g | ICIest=%7.4g ' ...
            '| maxU95(W/cluster)=%7.4g/%7.4g\n'], ...
            char(channels(c)),char(strategies(s)),meanBER,meanFER,meanSE,meanGoodput, ...
            meanPred,100*feas,meanPost,meanICI,maxW,maxB);
        row=table(channels(c),strategies(s),meanBER,meanFER,meanSE,meanGoodput, ...
            meanPred,feas,meanPost,meanICI,maxW,maxB, ...
            'VariableNames',{'Channel','Strategy','MeanBER','MeanFER','MeanNetSE_bpsHz', ...
            'MeanUncodedGoodput_bpsHz','MeanPredictedBER','PredictedFeasibleDataRatio', ...
            'MeanPosteriorVariance','MeanEstimatedICINoiseVariance', ...
            'MaxWilsonUpper95','MaxFrameBootstrapUpper95'});
        if isempty(summary), summary=row; else, summary=[summary;row]; end %#ok<AGROW>
    end
end

fprintf('\n--- B. Burst-adaptive mode and modulation mix by channel/SNR ---\n');
for c=1:numel(channels)
    snrs=unique(results.PilotSNR_dB(string(results.Channel)==channels(c)));
    for q=1:numel(snrs)
        mask=string(results.Channel)==channels(c) & string(results.Strategy)=='adaptive' & ...
            results.PilotSNR_dB==snrs(q);
        if ~any(mask), continue; end
        sm=mean(results.SpatialMultiplexingDataRatio(mask));
        al=mean(results.AlamoutiDataRatio(mask));
        qpsk=mean(results.QPSKDataRatio(mask));
        qam16=mean(results.QAM16DataRatio(mask));
        feasible=mean(results.PredictedTargetFeasibleDataRatio(mask));
        ber=mean(results.BER(mask));
        se=mean(results.NetSpectralEfficiency_bpsHz(mask));
        fprintf(['%-8s | SNR=%2g dB | SM=%4.0f%% Alamouti=%4.0f%% | QPSK=%4.0f%% 16QAM=%4.0f%% ' ...
            '| feasible=%4.0f%% | BER=%8.4g | SE=%5.2f\n'], ...
            char(channels(c)),snrs(q),100*sm,100*al,100*qpsk,100*qam16,100*feasible,ber,se);
    end
end

fprintf('\n--- C. Prediction calibration (predicted vs observed BER) ---\n');
for c=1:numel(channels)
    for s=1:numel(strategies)
        mask=string(results.Channel)==channels(c) & string(results.Strategy)==strategies(s);
        if ~any(mask), continue; end
        pred=mean(results.MeanPredictedBER(mask));
        obs=mean(results.BER(mask));
        fprintf('%-8s | %-22s | observed/predicted BER ratio = %.3g\n', ...
            char(channels(c)),char(strategies(s)),obs/max(pred,eps));
    end
end
fprintf(['\nNotes: BER is uncoded. Pilot-driven decisions are refreshed once per ' ...
    'pilot/data burst. Posterior channel-error variance and causal pilot-innovation ' ...
    'ICI estimates enter the MCS prediction and MMSE detector as equivalent ' ...
    'receive disturbance. Confidence maxima are worst-condition bounds, not pooled CIs.\n']);
end
