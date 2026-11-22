function plot_spatial_aware_lmmse_comparison(results,outDir)
%PLOT_SPATIAL_AWARE_LMMSE_COMPARISON Compare estimator BER and SE by profile.
channels = unique(string(results.Channel),'stable');
snrValues = unique(results.EsN0_dB,'sorted');
profiles = unique(string(results.CorrelationProfile),'stable');
estimators = unique(string(results.Estimator),'stable');
for c = 1:numel(channels)
    maskChannel = string(results.Channel)==channels(c);
    fig = figure('Visible','off','Name',sprintf('%s LMMSE comparison',char(channels(c))));
    subplot(2,1,1); hold on;
    for e = 1:numel(estimators)
        for s = 1:numel(snrValues)
            y = nan(numel(profiles),1);
            for p = 1:numel(profiles)
                idx = maskChannel & string(results.Estimator)==estimators(e) & ...
                    results.EsN0_dB==snrValues(s) & string(results.CorrelationProfile)==profiles(p);
                if any(idx), y(p)=results.BER(find(idx,1)); end
            end
            semilogy(1:numel(profiles),max(y,1e-8),'-o','LineWidth',1.1, ...
                'DisplayName',sprintf('%s, %g dB',char(estimators(e)),snrValues(s)));
        end
    end
    yline(1e-3,'k--','Target BER');
    grid on; ylabel('Measured BER');
    title(sprintf('%s: BER by estimator and spatial profile',upper(char(channels(c)))));
    legend('Location','eastoutside','Interpreter','none');
    subplot(2,1,2); hold on;
    for e = 1:numel(estimators)
        for s = 1:numel(snrValues)
            y = nan(numel(profiles),1);
            for p = 1:numel(profiles)
                idx = maskChannel & string(results.Estimator)==estimators(e) & ...
                    results.EsN0_dB==snrValues(s) & string(results.CorrelationProfile)==profiles(p);
                if any(idx), y(p)=results.SpectralEfficiency_bpsHz(find(idx,1)); end
            end
            plot(1:numel(profiles),y,'-o','LineWidth',1.1, ...
                'DisplayName',sprintf('%s, %g dB',char(estimators(e)),snrValues(s)));
        end
    end
    grid on; ylabel('Spectral efficiency (bit/s/Hz)'); xlabel('Correlation profile');
    xticks(1:numel(profiles)); xticklabels(profiles); xtickangle(20);
    legend('Location','eastoutside','Interpreter','none');
    saveas(fig,fullfile(outDir,sprintf('%s_spatial_aware_lmmse.png',char(channels(c)))));
    close(fig);
end
end
