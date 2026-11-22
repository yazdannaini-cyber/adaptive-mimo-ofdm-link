function plot_spatial_correlation_results(results,outDir)
%PLOT_SPATIAL_CORRELATION_RESULTS Plot BER and spectral efficiency vs correlation.
channels = unique(string(results.Channel),'stable');
snrValues = unique(results.EsN0_dB,'sorted');
profiles = unique(string(results.CorrelationProfile),'stable');
for c = 1:numel(channels)
    channelMask = string(results.Channel)==channels(c);
    fig = figure('Visible','off','Name',char(channels(c)+' correlation sweep'));
    subplot(2,1,1); hold on;
    for s = 1:numel(snrValues)
        y = nan(numel(profiles),1);
        for p = 1:numel(profiles)
            idx = channelMask & results.EsN0_dB==snrValues(s) & ...
                string(results.CorrelationProfile)==profiles(p);
            if any(idx), y(p)=results.BER(find(idx,1)); end
        end
        semilogy(1:numel(profiles),max(y,1e-7),'-o','LineWidth',1.2, ...
            'DisplayName',sprintf('%g dB',snrValues(s)));
    end
    yline(1e-3,'--','Target BER');
    grid on; ylabel('Measured BER'); title(sprintf('%s: BER vs spatial correlation',upper(char(channels(c)))));
    legend('Location','best');
    subplot(2,1,2); hold on;
    for s = 1:numel(snrValues)
        y = nan(numel(profiles),1);
        for p = 1:numel(profiles)
            idx = channelMask & results.EsN0_dB==snrValues(s) & ...
                string(results.CorrelationProfile)==profiles(p);
            if any(idx), y(p)=results.SpectralEfficiency_bpsHz(find(idx,1)); end
        end
        plot(1:numel(profiles),y,'-o','LineWidth',1.2, ...
            'DisplayName',sprintf('%g dB',snrValues(s)));
    end
    grid on; ylabel('bit/s/Hz'); xlabel('Correlation profile');
    xticks(1:numel(profiles)); xticklabels(profiles); xtickangle(20);
    legend('Location','best');
    saveas(fig,fullfile(outDir,sprintf('%s_spatial_correlation.png',char(channels(c)))));
    close(fig);
end
end
