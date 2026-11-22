function plot_pilot_mapped_adaptive_results(results,outDir)
%PLOT_PILOT_MAPPED_ADAPTIVE_RESULTS Create V15 BER and mode-mix plots.
if ~exist(outDir,'dir'), mkdir(outDir); end
channels = unique(string(results.Channel),'stable');
strategies = {'adaptive','fixed-sm-16qam','fixed-alamouti-qpsk'};
for c = 1:numel(channels)
    channelMask = string(results.Channel)==channels(c);
    snrs = unique(results.PilotSNR_dB(channelMask));
    fig = figure('Visible','off','Name',sprintf('%s V15 adaptive BER',char(channels(c))));
    hold on;
    for s = 1:numel(strategies)
        values = nan(size(snrs));
        for q = 1:numel(snrs)
            mask = channelMask & string(results.Strategy)==strategies{s} & results.PilotSNR_dB==snrs(q);
            if any(mask), values(q)=mean(results.BER(mask)); end
        end
        semilogy(snrs,max(values,1e-7),'-o','LineWidth',1.4,'DisplayName',strrep(strategies{s},'-',' '));
    end
    grid on; xlabel('Pilot/data SNR (dB)'); ylabel('Uncoded BER');
    title(sprintf('%s: adaptive versus fixed links',upper(char(channels(c)))),'Interpreter','none');
    legend('Location','best');
    saveas(fig,fullfile(outDir,sprintf('%s_v15_ber.png',char(channels(c)))));
    close(fig);

    fig = figure('Visible','off','Name',sprintf('%s V15 selected mode mix',char(channels(c))));
    sm = zeros(size(snrs)); al = zeros(size(snrs));
    for q = 1:numel(snrs)
        mask = channelMask & string(results.Strategy)=='adaptive' & results.PilotSNR_dB==snrs(q);
        if any(mask)
            sm(q)=mean(results.SpatialMultiplexingFrameRatio(mask));
            al(q)=mean(results.AlamoutiFrameRatio(mask));
        end
    end
    bar(snrs,[sm(:),al(:)],'stacked');
    ylim([0 1]); grid on; xlabel('Pilot/data SNR (dB)'); ylabel('Selected-frame fraction');
    legend({'Spatial Multiplexing','Alamouti STBC'},'Location','best');
    title(sprintf('%s: adaptive spatial mode selection',upper(char(channels(c)))),'Interpreter','none');
    saveas(fig,fullfile(outDir,sprintf('%s_v15_mode_mix.png',char(channels(c)))));
    close(fig);
end
end
