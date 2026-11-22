function summary = pilot_mapped_results_summary(results)
%PILOT_MAPPED_RESULTS_SUMMARY Print compact summaries of V13/V14 result rows.
% Usage after the V14 run:
%   summary = pilot_mapped_results_summary(results);
% Or with no argument, the saved CSV is read from pilot_mapped_results/.

if nargin < 1 || isempty(results)
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    csvPath = fullfile(projectRoot,'pilot_mapped_results','pilot_mapped_ofdm_validation.csv');
    if ~exist(csvPath,'file')
        error('pilot_mapped_results_summary:MissingCSV', ...
            'Run pilot_mapped_ofdm_validation_run first; CSV was not found.');
    end
    results = readtable(csvPath);
end
required = {'Channel','SpatialProfile','MaxDoppler_Hz','PilotSpacingSymbols', ...
    'Estimator','MeanNMSE','BER','FER','NetSpectralEfficiency_bpsHz','ICI_Ratio'};
if ~istable(results) || ~all(ismember(required,results.Properties.VariableNames))
    error('pilot_mapped_results_summary:InvalidTable', ...
        'Input table is missing one or more required pilot-mapped result columns.');
end

channels = unique(string(results.Channel),'stable');
estimators = unique(string(results.Estimator),'stable');
summary = struct();

fprintf('\n============================================================\n');
fprintf('PILOT-MAPPED OFDM RESULTS: %d rows\n',height(results));
fprintf('============================================================\n');
fprintf('\n--- A. Overall results by channel and estimator ---\n');
for c = 1:numel(channels)
    for m = 1:numel(estimators)
        mask = string(results.Channel)==channels(c) & string(results.Estimator)==estimators(m);
        part = results(mask,:);
        fprintf('%-8s | %-12s | NMSE=%8.4g | BER=%8.4g | FER=%6.3f | SE=%5.2f | ICI=%8.3g\n', ...
            char(channels(c)),char(estimators(m)),mean(part.MeanNMSE),mean(part.BER), ...
            mean(part.FER),mean(part.NetSpectralEfficiency_bpsHz),mean(part.ICI_Ratio));
    end
end

fprintf('\n--- B. Jakes-LMMSE NMSE and BER by Doppler (spacing 8, joint-high profile) ---\n');
for c = 1:numel(channels)
    for fd = unique(results.MaxDoppler_Hz).'
        mask = string(results.Channel)==channels(c) & ...
            string(results.SpatialProfile)=='joint-high' & ...
            results.PilotSpacingSymbols==8 & results.MaxDoppler_Hz==fd & ...
            string(results.Estimator)=='jakes-lmmse';
        part = results(mask,:);
        if ~isempty(part)
            fprintf('%-8s | fD=%3g Hz | NMSE=%8.4g | BER=%8.4g | FER=%6.3f | ICI=%8.3g\n', ...
                char(channels(c)),fd,mean(part.MeanNMSE),mean(part.BER), ...
                mean(part.FER),mean(part.ICI_Ratio));
        end
    end
end

fprintf('\n--- C. Jakes-LMMSE performance by pilot spacing ---\n');
for c = 1:numel(channels)
    for spacing = unique(results.PilotSpacingSymbols).'
        mask = string(results.Channel)==channels(c) & ...
            string(results.Estimator)=='jakes-lmmse' & ...
            results.PilotSpacingSymbols==spacing;
        part = results(mask,:);
        fprintf('%-8s | spacing=%2g | NMSE=%8.4g | BER=%8.4g | net SE=%5.2f | pilot overhead=%5.1f%%\n', ...
            char(channels(c)),spacing,mean(part.MeanNMSE),mean(part.BER), ...
            mean(part.NetSpectralEfficiency_bpsHz),100*mean(part.PilotOverheadRatio));
    end
end

fprintf('\n--- D. Measured ICI by channel, profile and Doppler (spacing 8) ---\n');
profileNames = unique(string(results.SpatialProfile),'stable');
for c = 1:numel(channels)
    for p = 1:numel(profileNames)
        for fd = unique(results.MaxDoppler_Hz).'
            mask = string(results.Channel)==channels(c) & ...
                string(results.SpatialProfile)==profileNames(p) & ...
                results.PilotSpacingSymbols==8 & results.MaxDoppler_Hz==fd & ...
                string(results.Estimator)=='perfect-csi';
            part = results(mask,:);
            if ~isempty(part)
                fprintf('%-8s | %-10s | fD=%3g Hz | ICI=%9.4g\n', ...
                    char(channels(c)),char(profileNames(p)),fd,mean(part.ICI_Ratio));
            end
        end
    end
end

fprintf('\n--- E. Diagonal perfect-CSI vs full-ICI oracle (spacing 8) ---\n');
for c = 1:numel(channels)
    for p = 1:numel(profileNames)
        for fd = unique(results.MaxDoppler_Hz).'
            common = string(results.Channel)==channels(c) & ...
                string(results.SpatialProfile)==profileNames(p) & ...
                results.PilotSpacingSymbols==8 & results.MaxDoppler_Hz==fd;
            diagPart = results(common & string(results.Estimator)=='perfect-csi',:);
            fullPart = results(common & string(results.Estimator)=='perfect-csi-ici',:);
            if ~isempty(diagPart) && ~isempty(fullPart)
                diagBER = mean(diagPart.BER);
                fullBER = mean(fullPart.BER);
                if diagBER > 0
                    reductionPct = 100*(diagBER-fullBER)/diagBER;
                else
                    reductionPct = NaN;
                end
                fprintf('%-8s | %-10s | fD=%3g Hz | diagonal BER=%8.4g | full-ICI BER=%8.4g | BER reduction=%7.2f%%\n', ...
                    char(channels(c)),char(profileNames(p)),fd,diagBER,fullBER,reductionPct);
            end
        end
    end
end

fprintf('\nV14 pilot-mapped results summary completed.\n');
end
