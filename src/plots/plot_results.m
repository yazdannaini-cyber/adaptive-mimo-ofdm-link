function plot_results(fixedTable, adaptiveTable, outputDir, saveFigures)
%PLOT_RESULTS Generate publication-style BER, NMSE, and adaptive-rate plots.
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
channels = unique(fixedTable.Channel);
for c = 1:numel(channels)
    channelName = channels(c);
    idx = fixedTable.Channel == channelName;
    T = fixedTable(idx,:);
    mods = unique(T.M);

    fig = figure('Visible','off','Color','w');
    hold on; grid on;
    lineStyles = {'-o','-s','-^','--d'};
    fields = {'BER_LS_MMSE','BER_LMMSE_MMSE','BER_PerfectCSI_MMSE','BER_LMMSE_ZF'};
    labels = {'LS + MMSE detector','LMMSE + MMSE detector','Perfect CSI + MMSE detector','LMMSE + ZF detector'};
    for m = 1:numel(mods)
        sub = T(T.M == mods(m),:);
        for q = 1:numel(fields)
            values = sub.(fields{q});
            values(values <= 0) = NaN; % avoid presenting zero-error estimates as exact zero BER
            semilogy(sub.EbN0_dB,values,lineStyles{q},'LineWidth',1.25, ...
                'DisplayName',sprintf('%d-QAM | %s',mods(m),labels{q}));
        end
    end
    xlabel('E_b/N_0 (dB)'); ylabel('Bit Error Rate (BER)');
    title(sprintf('Fixed-Modulation BER Performance | %s Channel',char(upper(channelName))));
    legend('Location','southwest'); ylim([1e-5 1]);
    fileName = fullfile(outputDir,sprintf('BER_Fixed_%s.png',char(channelName)));
    if saveFigures, print(fig,fileName,'-dpng','-r160'); end
    close(fig);

    % NMSE is independent of the selected data modulation; use the first M rows.
    sub = T(T.M == mods(1),:);
    fig = figure('Visible','off','Color','w');
    semilogy(sub.EbN0_dB,sub.NMSE_LS,'-o','LineWidth',1.5,'DisplayName','LS'); hold on;
    semilogy(sub.EbN0_dB,sub.NMSE_LMMSE,'-s','LineWidth',1.5,'DisplayName','LMMSE');
    grid on; xlabel('E_b/N_0 (dB)'); ylabel('Channel NMSE');
    title(sprintf('Channel Estimation Accuracy (4-QAM reference) | %s Channel',char(upper(channelName))));
    legend('Location','southwest');
    fileName = fullfile(outputDir,sprintf('NMSE_%s.png',char(channelName)));
    if saveFigures, print(fig,fileName,'-dpng','-r160'); end
    close(fig);
end

adaptiveChannels = unique(adaptiveTable.Channel);
for c = 1:numel(adaptiveChannels)
    channelName = adaptiveChannels(c);
    idx = adaptiveTable.Channel == channelName;
    T = adaptiveTable(idx,:);
    fig = figure('Visible','off','Color','w');
    yyaxis left;
    plot(T.EsN0_dB,T.SpectralEfficiency_bpsHz,'-o','LineWidth',1.5, ...
        'DisplayName','Spectral efficiency');
    ylabel('Spectral efficiency (bit/s/Hz)');
    yyaxis right;
    semilogy(T.EsN0_dB,max(T.BER,1e-7),'-s','LineWidth',1.5, ...
        'DisplayName','Measured BER'); hold on;
    if any(strcmp('PredictedBER_PayloadWeighted',T.Properties.VariableNames))
        semilogy(T.EsN0_dB,max(T.PredictedBER_PayloadWeighted,1e-7),'--^', ...
            'LineWidth',1.25,'DisplayName','Payload-weighted predicted BER');
    end
    if any(strcmp('TargetBER',T.Properties.VariableNames))
        yline(T.TargetBER(1),':','Target BER','LineWidth',1.25,'DisplayName','Target BER');
    end
    ylabel('BER');
    grid on; xlabel('E_s/N_0 (dB)');
    title(sprintf('Adaptive Modulation Performance | %s Channel',char(upper(channelName))));
    legend('Location','best');
    fileName = fullfile(outputDir,sprintf('Adaptive_Performance_%s.png',char(channelName)));
    if saveFigures, print(fig,fileName,'-dpng','-r160'); end
    close(fig);

    if any(strcmp('PredictedTargetMetFrameRatio',T.Properties.VariableNames))
        fig = figure('Visible','off','Color','w');
        plot(T.EsN0_dB,100*T.PredictedTargetMetFrameRatio,'-o','LineWidth',1.5, ...
            'DisplayName','Selected mode prediction meets target'); hold on;
        plot(T.EsN0_dB,100*T.FallbackQPSKFrameRatio,'--s','LineWidth',1.5, ...
            'DisplayName','QPSK fallback (target not guaranteed)');
        grid on; ylim([0 100]); xlabel('E_s/N_0 (dB)'); ylabel('Frames (%)');
        title(sprintf('Adaptive Link Feasibility | %s Channel',char(upper(channelName))));
        legend('Location','best');
        fileName = fullfile(outputDir,sprintf('Adaptive_Feasibility_%s.png',char(channelName)));
        if saveFigures, print(fig,fileName,'-dpng','-r160'); end
        close(fig);
    end
end
end
