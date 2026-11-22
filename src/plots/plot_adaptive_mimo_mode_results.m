function plot_adaptive_mimo_mode_results(results, outputDirectory)
%PLOT_ADAPTIVE_MIMO_MODE_RESULTS Plot joint mode/modulation adaptation results.
if ~istable(results) || isempty(results)
    error('plot_adaptive_mimo_mode_results:InvalidResults','results must be a nonempty table.');
end
if ~exist(outputDirectory,'dir'), mkdir(outputDirectory); end
channels = unique(string(results.Channel),'stable');

fig = figure('Visible','off'); hold on;
for c = 1:numel(channels)
    rows = results(string(results.Channel)==channels(c),:);
    [x,idx] = sort(rows.EsN0_dB);
    semilogy(x,max(rows.BER(idx),1e-7),'-o','LineWidth',1.5, ...
        'DisplayName',char(channels(c)));
end
yline(1e-3,'--','Target BER = 10^{-3}','DisplayName','Target BER');
grid on; xlabel('Es/N0 (dB)'); ylabel('Measured BER');
title('Joint MIMO-mode and modulation adaptation: BER');
legend('Location','southwest','Interpreter','none');
saveas(fig,fullfile(outputDirectory,'adaptive_mimo_mode_BER.png')); close(fig);

fig = figure('Visible','off'); hold on;
for c = 1:numel(channels)
    rows = results(string(results.Channel)==channels(c),:);
    [x,idx] = sort(rows.EsN0_dB);
    plot(x,rows.SpectralEfficiency_bpsHz(idx),'-o','LineWidth',1.5, ...
        'DisplayName',char(channels(c)));
end
grid on; xlabel('Es/N0 (dB)'); ylabel('Spectral efficiency (bit/s/Hz)');
title('Joint adaptation: net spectral efficiency');
legend('Location','best','Interpreter','none');
saveas(fig,fullfile(outputDirectory,'adaptive_mimo_mode_spectral_efficiency.png')); close(fig);

fig = figure('Visible','off'); hold on;
for c = 1:numel(channels)
    rows = results(string(results.Channel)==channels(c),:);
    [x,idx] = sort(rows.EsN0_dB);
    plot(x,rows.SpatialMultiplexingFrameRatio(idx),'-o','LineWidth',1.5, ...
        'DisplayName',[char(channels(c)) ' spatial']);
    plot(x,rows.AlamoutiFrameRatio(idx),'--s','LineWidth',1.5, ...
        'DisplayName',[char(channels(c)) ' Alamouti']);
end
grid on; xlabel('Es/N0 (dB)'); ylabel('Selected-frame fraction'); ylim([0 1]);
title('Adaptive MIMO transmission-mode selection');
legend('Location','best','Interpreter','none');
saveas(fig,fullfile(outputDirectory,'adaptive_mimo_mode_selection.png')); close(fig);
end
