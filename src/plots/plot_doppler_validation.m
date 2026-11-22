function plot_doppler_validation(results, acfResults, outDir)
%PLOT_DOPPLER_VALIDATION Plot ICI leakage and temporal ACF under Doppler.
if ~exist(outDir,'dir'), mkdir(outDir); end
channels = unique(string(results.Channel),'stable');
for c = 1:numel(channels)
    mask = string(results.Channel)==channels(c);
    fig = figure('Visible','off','Name',sprintf('%s Doppler validation',char(channels(c))));
    subplot(1,2,1);
    semilogy(results.MaxDoppler_Hz(mask),max(results.ICI_LeakageRatio(mask),1e-12),'-o','LineWidth',1.4);
    grid on; xlabel('Maximum Doppler (Hz)'); ylabel('ICI leakage ratio');
    title(sprintf('%s: ICI versus Doppler',upper(char(channels(c)))));
    subplot(1,2,2); hold on;
    fdValues = unique(acfResults.MaxDoppler_Hz(acfResults.Channel==channels(c)),'sorted');
    for k = 1:numel(fdValues)
        m = acfResults.Channel==channels(c) & acfResults.MaxDoppler_Hz==fdValues(k);
        plot(acfResults.LagSymbols(m),acfResults.EstimatedACF_Real(m),'-o', ...
            'DisplayName',sprintf('Estimated %g Hz',fdValues(k)));
        plot(acfResults.LagSymbols(m),acfResults.TheoreticalJakesACF(m),'--', ...
            'DisplayName',sprintf('Jakes %g Hz',fdValues(k)));
    end
    grid on; xlabel('Lag (OFDM symbols)'); ylabel('Normalized ACF');
    title('Diffuse fading temporal correlation');
    legend('Location','best');
    saveas(fig,fullfile(outDir,sprintf('%s_doppler_validation.png',char(channels(c)))));
    close(fig);
end
end
