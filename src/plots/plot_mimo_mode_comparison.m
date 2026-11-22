function plot_mimo_mode_comparison(results, outputDirectory)
%PLOT_MIMO_MODE_COMPARISON Save BER and spectral-efficiency comparison plots.
if ~istable(results) || isempty(results)
    error('plot_mimo_mode_comparison:InvalidResults','results must be a nonempty table.');
end
if ~exist(outputDirectory,'dir')
    mkdir(outputDirectory);
end
channels = unique(string(results.Channel),'stable');
modulations = unique(results.ModulationOrder,'stable');
for c = 1:numel(channels)
    for m = 1:numel(modulations)
        subset = results(string(results.Channel)==channels(c) & ...
            results.ModulationOrder==modulations(m),:);
        if isempty(subset)
            continue;
        end
        modes = unique(string(subset.TxMode),'stable');

        figBer = figure('Visible','off');
        hold on;
        for modeIndex = 1:numel(modes)
            modeRows = subset(string(subset.TxMode)==modes(modeIndex),:);
            [snrSorted,order] = sort(modeRows.EsN0_dB);
            berSorted = modeRows.BER(order);
            semilogy(snrSorted,max(berSorted,1e-7),'-o','LineWidth',1.4, ...
                'DisplayName',char(modes(modeIndex)));
        end
        yline(1e-3,'--','Target BER = 10^{-3}','DisplayName','Target BER');
        grid on;
        xlabel('Es/N0 (dB)');
        ylabel('Measured BER');
        title(sprintf('%s channel, %d-QAM: BER comparison',char(channels(c)),modulations(m)));
        legend('Location','southwest','Interpreter','none');
        saveas(figBer,fullfile(outputDirectory,sprintf('ber_%s_%dqam.png',char(channels(c)),modulations(m))));
        close(figBer);

        figFer = figure('Visible','off');
        hold on;
        for modeIndex = 1:numel(modes)
            modeRows = subset(string(subset.TxMode)==modes(modeIndex),:);
            [snrSorted,order] = sort(modeRows.EsN0_dB);
            ferSorted = modeRows.FER(order);
            semilogy(snrSorted,max(ferSorted,1e-7),'-o','LineWidth',1.4, ...
                'DisplayName',char(modes(modeIndex)));
        end
        grid on;
        xlabel('Es/N0 (dB)');
        ylabel('Frame error rate (FER)');
        title(sprintf('%s channel, %d-QAM: frame reliability',char(channels(c)),modulations(m)));
        legend('Location','southwest','Interpreter','none');
        saveas(figFer,fullfile(outputDirectory,sprintf('fer_%s_%dqam.png',char(channels(c)),modulations(m))));
        close(figFer);

        figSe = figure('Visible','off');
        hold on;
        for modeIndex = 1:numel(modes)
            modeRows = subset(string(subset.TxMode)==modes(modeIndex),:);
            [snrSorted,order] = sort(modeRows.EsN0_dB);
            seSorted = modeRows.SpectralEfficiency_bpsHz(order);
            plot(snrSorted,seSorted,'-o','LineWidth',1.4, ...
                'DisplayName',char(modes(modeIndex)));
        end
        grid on;
        xlabel('Es/N0 (dB)');
        ylabel('Spectral efficiency (bit/s/Hz)');
        title(sprintf('%s channel, %d-QAM: spectral efficiency',char(channels(c)),modulations(m)));
        legend('Location','best','Interpreter','none');
        saveas(figSe,fullfile(outputDirectory,sprintf('se_%s_%dqam.png',char(channels(c)),modulations(m))));
        close(figSe);
    end
end
end
