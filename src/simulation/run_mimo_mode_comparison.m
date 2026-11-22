function resultTable = run_mimo_mode_comparison(cfg, channelName)
%RUN_MIMO_MODE_COMPARISON Compare 2x2 spatial multiplexing and Alamouti STBC.
% Both modes use the same per-antenna symbol energy, pilot structure,
% channel realization, channel estimate, and per-subcarrier noise level.
% Spatial multiplexing carries two independent QAM streams per OFDM symbol.
% Alamouti STBC carries one QAM symbol per channel use (rate-1 code).
%
% The reported FER is the fraction of simulated data frames with >=1 bit error.
% No FEC decoder, CRC, or retransmission protocol is included in this comparison.

if ~isfield(cfg,'comparison')
    error('run_mimo_mode_comparison:MissingComparisonConfig', ...
        'cfg.comparison must define snrDb, modulations, and numFrames.');
end
if cfg.system.nTx ~= 2
    error('run_mimo_mode_comparison:TransmitAntennaCount', ...
        'This comparison currently requires exactly two transmit antennas.');
end
if mod(cfg.system.nDataSymbols,2) ~= 0
    error('run_mimo_mode_comparison:OddDataSymbols', ...
        'The number of data OFDM symbols must be even for Alamouti pairing.');
end

fprintf('\n[MIMO mode comparison] Channel: %s\n',upper(char(channelName)));
N = cfg.system.nFFT;
Nt = cfg.system.nTx;
Nr = cfg.system.nRx;
Ns = cfg.system.nDataSymbols;
cp = cfg.system.cpLength;
Ppilot = cfg.system.pilotPower;
numFrames = cfg.comparison.numFrames;
mods = cfg.comparison.modulations;
snrValues = cfg.comparison.snrDb;
[meanChannel,covarianceScale] = channel_estimation_prior(cfg,channelName);

if isempty(mods) || any(mods < 4) || any(abs(log2(mods)-round(log2(mods))) > eps) || ...
        any(abs(log2(sqrt(mods))-round(log2(sqrt(mods)))) > eps) || ...
        numFrames < 1 || numFrames ~= floor(numFrames)
    error('run_mimo_mode_comparison:InvalidSweep', 'Invalid modulation orders or frame count.');
end

resultTable = [];
for modulationIndex = 1:numel(mods)
    M = mods(modulationIndex);
    bitsPerQamSymbol = log2(M);
    for snrDb = snrValues
        noiseVarianceFD = 10^(-snrDb/10);
        noiseVarianceTime = noiseVarianceFD/N;
        channelErrorVariance = channel_estimation_error_variance( ...
            cfg.channel.pdp,noiseVarianceFD,Ppilot,N,covarianceScale);
        effectiveNoiseVariance = noiseVarianceFD + Nt*channelErrorVariance;

        totalBits = zeros(1,2);
        bitErrors = zeros(1,2);
        frameErrors = zeros(1,2);
        totalCorrectBits = zeros(1,2);
        frameBitErrors = cell(1,2);
        framePayloadBits = cell(1,2);
        for modeIndex = 1:2
            frameBitErrors{modeIndex} = zeros(numFrames,1);
            framePayloadBits{modeIndex} = zeros(numFrames,1);
        end

        for frameIndex = 1:numFrames
            taps = generate_mimo_channel(cfg,channelName);

            % Two orthogonal pilot OFDM symbols identify the two transmit links.
            Xpilot = complex(zeros(N,Nt,Nt));
            for txIndex = 1:Nt
                Xpilot(:,txIndex,txIndex) = sqrt(Ppilot);
            end
            pilotTx = ofdm_modulate(Xpilot,cp);
            pilotRx = apply_mimo_ofdm_channel(pilotTx,taps,noiseVarianceTime);
            pilotYfd = ofdm_demodulate(pilotRx,N,cp);
            Hls = channel_estimation_ls(pilotYfd,Ppilot);
            Hlmmse = channel_estimation_lmmse(Hls,cfg.channel.pdp,noiseVarianceFD, ...
                Ppilot,meanChannel,covarianceScale);

            % Mode 1: 2x2 spatial multiplexing with linear MMSE detection.
            bitsSpatial = randi([0 1],N*Nt*Ns*bitsPerQamSymbol,1);
            Xspatial = reshape(qam_mod(bitsSpatial,M),N,Nt,Ns);
            txSpatial = ofdm_modulate(Xspatial,cp);
            rxSpatial = apply_mimo_ofdm_channel(txSpatial,taps,noiseVarianceTime);
            Yspatial = ofdm_demodulate(rxSpatial,N,cp);
            [XspatialHat,~] = mimo_equalize(Yspatial,Hlmmse,effectiveNoiseVariance,'mmse');
            bitsSpatialHat = qam_demod(XspatialHat(:),M);
            errorsSpatial = sum(bitsSpatialHat ~= bitsSpatial);
            nBitsSpatial = numel(bitsSpatial);

            totalBits(1) = totalBits(1)+nBitsSpatial;
            bitErrors(1) = bitErrors(1)+errorsSpatial;
            frameErrors(1) = frameErrors(1)+double(errorsSpatial > 0);
            totalCorrectBits(1) = totalCorrectBits(1)+(nBitsSpatial-errorsSpatial);
            frameBitErrors{1}(frameIndex) = errorsSpatial;
            framePayloadBits{1}(frameIndex) = nBitsSpatial;

            % Mode 2: rate-1 Alamouti STBC over adjacent OFDM symbols.
            bitsStbc = randi([0 1],N*Ns*bitsPerQamSymbol,1);
            Xinformation = reshape(qam_mod(bitsStbc,M),N,2,Ns/2);
            Xstbc = alamouti_encode(Xinformation);
            txStbc = ofdm_modulate(Xstbc,cp);
            rxStbc = apply_mimo_ofdm_channel(txStbc,taps,noiseVarianceTime);
            Ystbc = ofdm_demodulate(rxStbc,N,cp);
            [XstbcHat,~] = alamouti_decode(Ystbc,Hlmmse,effectiveNoiseVariance);
            bitsStbcHat = qam_demod(XstbcHat(:),M);
            errorsStbc = sum(bitsStbcHat ~= bitsStbc);
            nBitsStbc = numel(bitsStbc);

            totalBits(2) = totalBits(2)+nBitsStbc;
            bitErrors(2) = bitErrors(2)+errorsStbc;
            frameErrors(2) = frameErrors(2)+double(errorsStbc > 0);
            totalCorrectBits(2) = totalCorrectBits(2)+(nBitsStbc-errorsStbc);
            frameBitErrors{2}(frameIndex) = errorsStbc;
            framePayloadBits{2}(frameIndex) = nBitsStbc;
        end

        modeNames = {'spatial-multiplexing','alamouti-stbc'};
        for modeIndex = 1:2
            metrics = bit_error_metrics(bitErrors(modeIndex),totalBits(modeIndex), ...
                frameErrors(modeIndex),numFrames);
            ferMetrics = bit_error_metrics(frameErrors(modeIndex),numFrames,0,1);
            bootstrapSeed = cfg.statistics.bootstrapSeed + round(snrDb*100) + ...
                M*100 + modeIndex*10000 + sum(double(char(channelName)));
            clusterCI = cluster_bootstrap_ber_ci(frameBitErrors{modeIndex}, ...
                framePayloadBits{modeIndex},cfg.statistics.bootstrapSamples, ...
                cfg.statistics.confidenceLevel,bootstrapSeed);

            if modeIndex == 1
                streamCount = Nt;
            else
                streamCount = 1; % Alamouti rate is one QAM symbol per channel use.
            end
            spectralEfficiency = streamCount*bitsPerQamSymbol*Ns/(Ns+Nt)*N/(N+cp);
            goodput = (totalCorrectBits(modeIndex)/numFrames)/((Ns+Nt)*N)*N/(N+cp);

            row = table(string(channelName),snrDb,M,string(modeNames{modeIndex}), ...
                metrics.BER,metrics.BLER,metrics.BER_CI95_Lower,metrics.BER_CI95_Upper, ...
                ferMetrics.BER_CI95_Lower,ferMetrics.BER_CI95_Upper,frameErrors(modeIndex), ...
                clusterCI(1),clusterCI(2),spectralEfficiency,goodput, ...
                totalBits(modeIndex),numFrames, ...
                'VariableNames',{'Channel','EsN0_dB','ModulationOrder','TxMode', ...
                'BER','FER','BER_Wilson_CI95_Lower','BER_Wilson_CI95_Upper', ...
                'FER_Wilson_CI95_Lower','FER_Wilson_CI95_Upper','FrameErrorFrames', ...
                'BER_FrameBootstrap_CI95_Lower','BER_FrameBootstrap_CI95_Upper', ...
                'SpectralEfficiency_bpsHz','Goodput_bpsHz','PayloadBits','NumFrames'});
            if isempty(resultTable)
                resultTable = row;
            else
                resultTable = [resultTable;row]; %#ok<AGROW>
            end

            fprintf(['  M=%2d | Es/N0=%2d dB | %-20s | BER=%.3g | FER=%.3g [%.2g, %.2g] | ' ...
                'SE=%.2f bit/s/Hz | Goodput=%.2f bit/s/Hz\n'], ...
                M,snrDb,modeNames{modeIndex},metrics.BER,metrics.BLER, ...
                ferMetrics.BER_CI95_Lower,ferMetrics.BER_CI95_Upper,spectralEfficiency,goodput);
        end
    end
end
end
