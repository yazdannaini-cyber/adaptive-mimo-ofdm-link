function resultTable = run_fixed_modulation(cfg, channelName)
%RUN_FIXED_MODULATION Compare channel estimation and MIMO detectors.
% BER is evaluated versus Eb/N0 for fixed QAM orders.

fprintf('\n[Fixed modulation] Channel: %s\n', upper(channelName));
resultTable = [];
N = cfg.system.nFFT;
Nt = cfg.system.nTx;
Nr = cfg.system.nRx;
Ns = cfg.system.nDataSymbols;
cp = cfg.system.cpLength;
Ppilot = cfg.system.pilotPower;
numFrames = cfg.simulation.numFrames;
[meanChannel,covarianceScale] = channel_estimation_prior(cfg,channelName);

for M = cfg.simulation.modulations
    bitsPerSymbol = log2(M);
    for ebnoDb = cfg.simulation.ebn0dB
        ebnoLinear = 10^(ebnoDb/10);
        noiseVarianceFD = 1/(bitsPerSymbol*ebnoLinear);
        noiseVarianceTime = noiseVarianceFD/N;
        % Variants: LS+MMSE, LMMSE+MMSE, perfect-CSI reference, LMMSE+ZF.
        bitErrors = zeros(1,4);
        blockErrors = zeros(1,4);
        frameBitErrors = zeros(numFrames,4);
        framePayloadBits = zeros(numFrames,1);
        totalBits = 0;
        nmseLS_num = 0;
        nmseLMMSE_num = 0;
        nmseDen = 0;

        for frame = 1:numFrames
            taps = generate_mimo_channel(cfg, channelName);
            Htrue = fft(taps, N, 1);

            % Orthogonal pilot OFDM symbols: one active Tx antenna per symbol.
            Xpilot = complex(zeros(N,Nt,Nt));
            for tx = 1:Nt
                Xpilot(:,tx,tx) = sqrt(Ppilot);
            end
            pilotTx = ofdm_modulate(Xpilot, cp);
            pilotRx = apply_mimo_ofdm_channel(pilotTx, taps, noiseVarianceTime);
            pilotYfd = ofdm_demodulate(pilotRx, N, cp);
            Hls = channel_estimation_ls(pilotYfd, Ppilot);
            Hlmmse = channel_estimation_lmmse(Hls,cfg.channel.pdp,noiseVarianceFD,Ppilot, ...
                meanChannel,covarianceScale);

            % Generate one common payload for a fair detector/estimator comparison.
            txBits = randi([0 1], N*Nt*Ns*bitsPerSymbol, 1);
            txSymbols = qam_mod(txBits, M);
            Xdata = reshape(txSymbols, N, Nt, Ns);
            dataTx = ofdm_modulate(Xdata, cp);
            dataRx = apply_mimo_ofdm_channel(dataTx, taps, noiseVarianceTime);
            Ydata = ofdm_demodulate(dataRx, N, cp);

            Hset = {Hls, Hlmmse, Htrue, Hlmmse};
            methods = {'mmse','mmse','mmse','zf'};
            for variant = 1:4
                [Xhat, ~] = mimo_equalize(Ydata, Hset{variant}, noiseVarianceFD, methods{variant});
                rxBits = qam_demod(Xhat(:), M);
                errors = sum(rxBits ~= txBits);
                bitErrors(variant) = bitErrors(variant) + errors;
                frameBitErrors(frame,variant) = errors;
                blockErrors(variant) = blockErrors(variant) + double(errors > 0);
            end
            framePayloadBits(frame) = numel(txBits);
            totalBits = totalBits + numel(txBits);
            nmseLS_num = nmseLS_num + sum(abs(Hls(:)-Htrue(:)).^2);
            nmseLMMSE_num = nmseLMMSE_num + sum(abs(Hlmmse(:)-Htrue(:)).^2);
            nmseDen = nmseDen + sum(abs(Htrue(:)).^2);
        end

        m1 = bit_error_metrics(bitErrors(1), totalBits, blockErrors(1), numFrames);
        m2 = bit_error_metrics(bitErrors(2), totalBits, blockErrors(2), numFrames);
        m3 = bit_error_metrics(bitErrors(3), totalBits, blockErrors(3), numFrames);
        m4 = bit_error_metrics(bitErrors(4), totalBits, blockErrors(4), numFrames);
        clusterCI = zeros(4,2);
        for variant = 1:4
            clusterCI(variant,:) = cluster_bootstrap_ber_ci(frameBitErrors(:,variant), ...
                framePayloadBits,cfg.statistics.bootstrapSamples,cfg.statistics.confidenceLevel, ...
                cfg.statistics.bootstrapSeed + M*1000 + ebnoDb*10 + variant + sum(double(channelName)));
        end
        row = table(string(channelName), M, string(sprintf('%d-QAM',M)), ebnoDb, ...
            m1.BER, m2.BER, m3.BER, m4.BER, ...
            m1.BLER, m2.BLER, m3.BLER, m4.BLER, ...
            m1.BER_CI95_Lower, m1.BER_CI95_Upper, ...
            m2.BER_CI95_Lower, m2.BER_CI95_Upper, ...
            clusterCI(1,1),clusterCI(2,1),clusterCI(3,1),clusterCI(4,1), ...
            clusterCI(1,2),clusterCI(2,2),clusterCI(3,2),clusterCI(4,2), ...
            nmseLS_num/max(nmseDen,eps), nmseLMMSE_num/max(nmseDen,eps), ...
            totalBits, numFrames, ...
            'VariableNames', {'Channel','M','Modulation','EbN0_dB', ...
            'BER_LS_MMSE','BER_LMMSE_MMSE','BER_PerfectCSI_MMSE','BER_LMMSE_ZF', ...
            'BLER_LS_MMSE','BLER_LMMSE_MMSE','BLER_PerfectCSI_MMSE','BLER_LMMSE_ZF', ...
            'BER_LS_CI95_Lower','BER_LS_CI95_Upper','BER_LMMSE_CI95_Lower','BER_LMMSE_CI95_Upper', ...
            'BER_LS_MMSE_ClusterCI95_Lower','BER_LMMSE_MMSE_ClusterCI95_Lower', ...
            'BER_PerfectCSI_MMSE_ClusterCI95_Lower','BER_LMMSE_ZF_ClusterCI95_Lower', ...
            'BER_LS_MMSE_ClusterCI95_Upper','BER_LMMSE_MMSE_ClusterCI95_Upper', ...
            'BER_PerfectCSI_MMSE_ClusterCI95_Upper','BER_LMMSE_ZF_ClusterCI95_Upper', ...
            'NMSE_LS','NMSE_LMMSE','TotalBits','NumFrames'});
        if isempty(resultTable)
            resultTable = row;
        else
            resultTable = [resultTable; row]; %#ok<AGROW>
        end
        fprintf('  %2d-QAM | Eb/N0=%2d dB | BER (LMMSE+MMSE)=%.3g | NMSE=%.3g\n', ...
            M, ebnoDb, m2.BER, nmseLMMSE_num/max(nmseDen,eps));
    end
end
end
