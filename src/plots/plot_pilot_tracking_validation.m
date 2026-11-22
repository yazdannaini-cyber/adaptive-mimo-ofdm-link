function plot_pilot_tracking_validation(results,outDir)
%PLOT_PILOT_TRACKING_VALIDATION Plot Doppler and pilot-spacing tracking studies.
if ~exist(outDir,'dir'), mkdir(outDir); end
channels = unique(string(results.Channel),'stable');
profiles = unique(string(results.SpatialProfile),'stable');
methods = unique(string(results.Estimator),'stable');

for c = 1:numel(channels)
    fig = figure('Visible','off','Name',sprintf('%s pilot tracking versus Doppler',char(channels(c))));
    for p = 1:numel(profiles)
        subplot(1,numel(profiles),p); hold on;
        for m = 1:numel(methods)
            mask = string(results.Channel)==channels(c) & ...
                string(results.SpatialProfile)==profiles(p) & ...
                results.PilotSpacingSymbols==8 & results.PilotSNR_dB==10 & ...
                string(results.Estimator)==methods(m);
            part = sortrows(results(mask,:),'MaxDoppler_Hz');
            if ~isempty(part)
                semilogy(part.MaxDoppler_Hz,max(part.MeanNMSE,1e-12),'-o', ...
                    'LineWidth',1.3,'DisplayName',char(methods(m)));
            end
        end
        grid on; xlabel('Maximum Doppler (Hz)'); ylabel('Data-symbol NMSE');
        title(char(profiles(p)),'Interpreter','none');
        if p==1, legend('Location','best'); end
    end
    sgtitle(sprintf('%s: tracking NMSE vs Doppler (pilot spacing 8, pilot SNR 10 dB)', ...
        upper(char(channels(c)))));
    saveas(fig,fullfile(outDir,sprintf('%s_pilot_tracking_vs_doppler.png',char(channels(c)))));
    close(fig);

    fig = figure('Visible','off','Name',sprintf('%s pilot tracking versus spacing',char(channels(c))));
    for p = 1:numel(profiles)
        subplot(1,numel(profiles),p); hold on;
        for m = 1:numel(methods)
            mask = string(results.Channel)==channels(c) & ...
                string(results.SpatialProfile)==profiles(p) & ...
                results.MaxDoppler_Hz==120 & results.PilotSNR_dB==10 & ...
                string(results.Estimator)==methods(m);
            part = sortrows(results(mask,:),'PilotSpacingSymbols');
            if ~isempty(part)
                semilogy(part.PilotSpacingSymbols,max(part.MeanNMSE,1e-12),'-o', ...
                    'LineWidth',1.3,'DisplayName',char(methods(m)));
            end
        end
        grid on; xlabel('Pilot update period (OFDM symbols)'); ylabel('Data-symbol NMSE');
        title(char(profiles(p)),'Interpreter','none');
        if p==1, legend('Location','best'); end
    end
    sgtitle(sprintf('%s: tracking NMSE vs pilot spacing (Doppler 120 Hz, pilot SNR 10 dB)', ...
        upper(char(channels(c)))));
    saveas(fig,fullfile(outDir,sprintf('%s_pilot_tracking_vs_spacing.png',char(channels(c)))));
    close(fig);
end
end
