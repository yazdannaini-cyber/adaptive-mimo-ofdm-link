function plot_pilot_mapped_ofdm_results(results,outDir)
%PLOT_PILOT_MAPPED_OFDM_RESULTS Plot actual-pilot OFDM tracking and ICI results.
if ~exist(outDir,'dir'), mkdir(outDir); end
channels = unique(string(results.Channel),'stable');
profiles = unique(string(results.SpatialProfile),'stable');
methods = unique(string(results.Estimator),'stable');

for c = 1:numel(channels)
    fig = figure('Visible','off','Name',sprintf('%s mapped-pilot NMSE',char(channels(c))));
    for p = 1:numel(profiles)
        subplot(1,numel(profiles),p); hold on;
        for m = 1:numel(methods)
            mask = string(results.Channel)==channels(c) & ...
                string(results.SpatialProfile)==profiles(p) & ...
                results.PilotSpacingSymbols==8 & string(results.Estimator)==methods(m);
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
    sgtitle(sprintf('%s: mapped-pilot channel tracking (spacing 8)',upper(char(channels(c)))));
    saveas(fig,fullfile(outDir,sprintf('%s_mapped_pilot_nmse.png',char(channels(c)))));
    close(fig);

    fig = figure('Visible','off','Name',sprintf('%s mapped-pilot BER',char(channels(c))));
    for p = 1:numel(profiles)
        subplot(1,numel(profiles),p); hold on;
        for m = 1:numel(methods)
            mask = string(results.Channel)==channels(c) & ...
                string(results.SpatialProfile)==profiles(p) & ...
                results.PilotSpacingSymbols==8 & string(results.Estimator)==methods(m);
            part = sortrows(results(mask,:),'MaxDoppler_Hz');
            if ~isempty(part)
                semilogy(part.MaxDoppler_Hz,max(part.BER,1e-7),'-o', ...
                    'LineWidth',1.3,'DisplayName',char(methods(m)));
            end
        end
        yline(1e-3,'--','Target BER 10^{-3}');
        grid on; xlabel('Maximum Doppler (Hz)'); ylabel('Measured BER');
        title(char(profiles(p)),'Interpreter','none');
        if p==1, legend('Location','best'); end
    end
    sgtitle(sprintf('%s: uncoded BER with mapped pilots (spacing 8)',upper(char(channels(c)))));
    saveas(fig,fullfile(outDir,sprintf('%s_mapped_pilot_ber.png',char(channels(c)))));
    close(fig);

    mask = string(results.Channel)==channels(c) & string(results.Estimator)=='perfect-csi' & ...
        results.PilotSpacingSymbols==8;
    iciRows = unique(results(mask,{'SpatialProfile','MaxDoppler_Hz','ICI_Ratio'}));
    fig = figure('Visible','off','Name',sprintf('%s mapped-pilot ICI',char(channels(c))));
    hold on;
    for p = 1:numel(profiles)
        part = sortrows(iciRows(string(iciRows.SpatialProfile)==profiles(p),:),'MaxDoppler_Hz');
        if ~isempty(part)
            plot(part.MaxDoppler_Hz,part.ICI_Ratio,'-o', ...
                'LineWidth',1.3,'DisplayName',char(profiles(p)));
        end
    end
    grid on; xlabel('Maximum Doppler (Hz)'); ylabel('Measured ICI / diagonal-channel power');
    legend('Location','best');
    title(sprintf('%s: ICI caused by within-symbol channel variation',upper(char(channels(c)))));
    saveas(fig,fullfile(outDir,sprintf('%s_mapped_pilot_ici.png',char(channels(c)))));
    close(fig);
end
end
