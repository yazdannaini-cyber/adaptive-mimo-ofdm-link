function plot_channel_estimator_mse_validation(results,outDir)
%PLOT_CHANNEL_ESTIMATOR_MSE_VALIDATION Plot direct estimator NMSE comparison.
channels = unique(string(results.Channel),'stable');
profiles = unique(string(results.CorrelationProfile),'stable');
for c = 1:numel(channels)
    figure('Visible','off','Name',char(channels(c)));
    for p = 1:numel(profiles)
        subset = results(results.Channel==channels(c) & ...
            results.CorrelationProfile==profiles(p),:);
        [snrValues,order] = sort(subset.PilotSNR_dB);
        perLink = subset.PerLink_NMSE(order);
        joint = subset.Joint_NMSE(order);
        subplot(2,3,p);
        semilogy(snrValues,perLink,'-o','LineWidth',1.2); hold on;
        semilogy(snrValues,joint,'-s','LineWidth',1.2);
        grid on;
        xlabel('Pilot SNR (dB)'); ylabel('NMSE');
        title(char(profiles(p)),'Interpreter','none');
        if p==1, legend('Per-link LMMSE','Joint LMMSE','Location','southwest'); end
    end
    sgtitle(sprintf('%s: Direct Channel-Estimation NMSE',upper(char(channels(c)))));
    saveas(gcf,fullfile(outDir,sprintf('%s_estimator_nmse.png',char(channels(c)))));
    close(gcf);
end

figure('Visible','off','Name','Joint estimator improvement');
for c = 1:numel(channels)
    subplot(1,numel(channels),c);
    subset = results(results.Channel==channels(c),:);
    groupedProfiles = unique(string(subset.CorrelationProfile),'stable');
    hold on;
    for p = 1:numel(groupedProfiles)
        q = subset(subset.CorrelationProfile==groupedProfiles(p),:);
        [snrValues,order] = sort(q.PilotSNR_dB);
        plot(snrValues,q.Joint_NMSE_Improvement_Percent(order),'-o','LineWidth',1.1);
    end
    yline(0,'--'); grid on;
    xlabel('Pilot SNR (dB)'); ylabel('NMSE improvement (%)');
    title(upper(char(channels(c))));
    legend(cellstr(groupedProfiles),'Interpreter','none','Location','best');
end
saveas(gcf,fullfile(outDir,'joint_lmmse_nmse_improvement.png'));
close(gcf);
end
