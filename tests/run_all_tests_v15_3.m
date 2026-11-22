function run_all_tests_v15_3()
%RUN_ALL_TESTS_V15_3 Regression tests plus joint MIMO/MCS/FEC controller.

rng(1234,'twister');
fprintf('\nRunning MIMO-OFDM project tests...\n');

% Test 1: exact hard-decision QAM round trip in the absence of noise.
for M = [4 16 64]
    bits = randi([0 1], 600*log2(M), 1);
    recovered = qam_demod(qam_mod(bits,M),M);
    assert(isequal(bits,recovered),sprintf('QAM round-trip failed for M=%d.',M));
end
fprintf('  Test 1 passed: Gray-coded QAM mapping/demapping.\n');

% Test 2: OFDM modulation/demodulation round trip.
N = 64; Nt = 2; Ns = 3; cp = 16;
X = randn(N,Nt,Ns) + 1i*randn(N,Nt,Ns);
tx = ofdm_modulate(X,cp);
Y = ofdm_demodulate(tx,N,cp);
assert(max(abs(X(:)-Y(:))) < 1e-10,'OFDM round-trip failed.');
fprintf('  Test 2 passed: IFFT/CP/FFT round trip.\n');

% Test 3: noiseless AWGN-reference MIMO channel behaves as identity.
cfg = default_config();
cfg.system.nFFT = N; cfg.system.cpLength = cp;
cfg.system.nTx = Nt; cfg.system.nRx = Nt;
cfg.channel.numTaps = 1; cfg.channel.pdp = 1;
validate_config(cfg);
taps = generate_mimo_channel(cfg,'awgn');
rx = apply_mimo_ofdm_channel(tx,taps,0);
Ych = ofdm_demodulate(rx,N,cp);
assert(max(abs(X(:)-Ych(:))) < 1e-10,'Noiseless AWGN-reference channel test failed.');
fprintf('  Test 3 passed: OFDM signal through identity channel.\n');

% Test 4: multipath channel with sufficient CP matches the per-subcarrier model.
cfgRay = default_config();
tapsRay = generate_mimo_channel(cfgRay,'rayleigh');
Xray = randn(N,Nt,Ns) + 1i*randn(N,Nt,Ns);
txRay = ofdm_modulate(Xray,cp);
rxRay = apply_mimo_ofdm_channel(txRay,tapsRay,0);
Yray = ofdm_demodulate(rxRay,N,cp);
Hray = fft(tapsRay,N,1);
Yexpected = complex(zeros(N,Nt,Ns));
for ss = 1:Ns
    for kk = 1:N
        xk = reshape(Xray(kk,:,ss),Nt,1);
        for rr = 1:Nt
            hrow = reshape(Hray(kk,rr,:),1,Nt);
            Yexpected(kk,rr,ss) = hrow*xk;
        end
    end
end
assert(max(abs(Yray(:)-Yexpected(:))) < 1e-10,'Multipath OFDM frequency-domain relation failed.');
fprintf('  Test 4 passed: multipath channel and cyclic-prefix circular-convolution property.\n');

% Test 5: LS pilot estimate equals the true frequency response without noise.
Xpilot = complex(zeros(N,Nt,Nt));
for t = 1:Nt, Xpilot(:,t,t) = 1; end
pilotTx = ofdm_modulate(Xpilot,cp);
pilotRx = apply_mimo_ofdm_channel(pilotTx,taps,0);
pilotY = ofdm_demodulate(pilotRx,N,cp);
Hls = channel_estimation_ls(pilotY,1);
Htrue = fft(taps,N,1);
assert(max(abs(Hls(:)-Htrue(:))) < 1e-10,'Noiseless LS channel estimate failed.');
fprintf('  Test 5 passed: LS channel estimation with orthogonal pilots.\n');

% Test 6: frequency-domain LMMSE estimate is exact for a noiseless valid channel.
Hlm = channel_estimation_lmmse(Hls,1,0,1);
assert(max(abs(Hlm(:)-Htrue(:))) < 1e-8,'Noiseless LMMSE channel estimate failed.');
fprintf('  Test 6 passed: LMMSE smoothing in the noiseless case.\n');

% Test 7: ZF and low-noise MMSE detector recover a well-conditioned MIMO stream.
H = [1.0 0.25+0.1i; 0.15-0.1i 0.9];
X0 = randn(N,Nt,1) + 1i*randn(N,Nt,1);
Harray = repmat(reshape(H,1,Nt,Nt),N,1,1);
Y0 = complex(zeros(N,Nt,1));
for k = 1:N, Y0(k,:,1) = (H*reshape(X0(k,:,1),Nt,1)).'; end
[Xzf,~] = mimo_equalize(Y0,Harray,1e-12,'zf');
[Xmmse,~] = mimo_equalize(Y0,Harray,1e-12,'mmse');
assert(max(abs(Xzf(:)-X0(:))) < 1e-8,'ZF detector recovery test failed.');
assert(max(abs(Xmmse(:)-X0(:))) < 1e-7,'MMSE detector recovery test failed.');
fprintf('  Test 7 passed: ZF/MMSE MIMO detection.\n');

% Test 8: metrics return correct BER/BLER and sensible confidence bounds.
metric = bit_error_metrics(10,1000,2,20);
assert(abs(metric.BER-0.01) < eps && abs(metric.BLER-0.1) < eps, 'Error metric calculation failed.');
assert(metric.BER_CI95_Lower <= metric.BER && metric.BER <= metric.BER_CI95_Upper, ...
    'Wilson interval does not contain the measured BER.');
fprintf('  Test 8 passed: BER/BLER and Wilson confidence interval.\n');

% Test 9: model-aware Rician LMMSE recovers a valid noiseless channel.
cfgR = default_config();
[meanRician,covScaleRician] = channel_estimation_prior(cfgR,'rician');
tapsR = generate_mimo_channel(cfgR,'rician');
HtrueR = fft(tapsR,N,1);
HricianEstimate = channel_estimation_lmmse(HtrueR,cfgR.channel.pdp,0, ...
    cfgR.system.pilotPower,meanRician,covScaleRician);
assert(max(abs(HricianEstimate(:)-HtrueR(:))) < 1e-8, ...
    'Mean-aware noiseless Rician LMMSE estimate failed.');
fprintf('  Test 9 passed: mean-aware LMMSE with Rician LOS prior.\n');

% Test 10: posterior channel-estimation variance decreases with lower noise
% and with the weaker diffuse component of a Rician channel.
varHighNoise = channel_estimation_error_variance(cfgR.channel.pdp,0.05, ...
    cfgR.system.pilotPower,N,1);
varLowNoise = channel_estimation_error_variance(cfgR.channel.pdp,0.005, ...
    cfgR.system.pilotPower,N,1);
varRician = channel_estimation_error_variance(cfgR.channel.pdp,0.05, ...
    cfgR.system.pilotPower,N,covScaleRician);
assert(all(varHighNoise >= 0) && mean(varLowNoise) < mean(varHighNoise), ...
    'Posterior channel variance must decrease as pilot noise decreases.');
assert(mean(varRician) < mean(varHighNoise), ...
    'Rician diffuse channel uncertainty should be below Rayleigh uncertainty.');
fprintf('  Test 10 passed: posterior channel-estimation uncertainty.\n');

% Test 11: vector noise input is equivalent to scalar noise when constant.
noiseScalar = 0.03;
[Xscalar,Sscalar] = mimo_equalize(Y0,Harray,noiseScalar,'mmse');
[Xvector,Svector] = mimo_equalize(Y0,Harray,noiseScalar*ones(N,1),'mmse');
assert(max(abs(Xscalar(:)-Xvector(:))) < 1e-12 && ...
    max(abs(Sscalar(:)-Svector(:))) < 1e-12, ...
    'Per-subcarrier noise vector handling is inconsistent.');
fprintf('  Test 11 passed: per-subcarrier effective-noise support.\n');

% Test 12: the QAM BER predictor is monotonic and penalizes higher QAM order.
berQpskLow = qam_ber_approx(5,4);
berQpskHigh = qam_ber_approx(20,4);
ber16AtSameSinr = qam_ber_approx(20,16);
assert(berQpskHigh < berQpskLow && ber16AtSameSinr > berQpskHigh, ...
    'QAM BER predictor monotonicity/modulation-order test failed.');
fprintf('  Test 12 passed: post-SINR QAM BER prediction.\n');

% Test 13: the AWGN reference prior matches its deterministic identity channel.
[meanAwgn,covScaleAwgn] = channel_estimation_prior(cfg,'awgn');
Hidentity = fft(taps,N,1);
assert(covScaleAwgn == 0 && max(abs(meanAwgn(:)-Hidentity(:))) < 1e-12, ...
    'AWGN channel prior must match the deterministic identity response.');
HnoisyDummy = Hidentity+0.1*(randn(size(Hidentity))+1i*randn(size(Hidentity)));
HawgnEstimate = channel_estimation_lmmse(HnoisyDummy,cfg.channel.pdp,0.1, ...
    cfg.system.pilotPower,meanAwgn,covScaleAwgn);
assert(max(abs(HawgnEstimate(:)-Hidentity(:))) < 1e-12, ...
    'Deterministic AWGN prior must return the known identity response.');
fprintf('  Test 13 passed: deterministic AWGN matched prior.\n');


% Test 14: predicted-BER selector chooses highest eligible mode and clearly
% reports when QPSK is only a target-violating fallback.
[idxEligible,feasibleEligible] = select_qam_by_predicted_ber([2e-2 8e-4 2e-4],1e-3);
[idxFallback,feasibleFallback] = select_qam_by_predicted_ber([4e-2 2e-2 1e-2],1e-3);
assert(idxEligible == 3 && feasibleEligible, ...
    'Adaptive selector did not choose the highest target-compliant mode.');
assert(idxFallback == 1 && ~feasibleFallback, ...
    'Adaptive selector must flag QPSK fallback as target-infeasible.');
fprintf('  Test 14 passed: adaptive mode selection and fallback semantics.\n');

% Test 15: aggregate predicted BER is weighted by each frame's payload bits;
% frame feasibility and measured target compliance remain distinct metrics.
stats = adaptive_link_metrics([1e-2;1e-3],[100;900],[false;true],[true;false], ...
    1.5e-3,2.1e-3,1e-3);
assert(abs(stats.PredictedBER_PayloadWeighted-1.9e-3) < 1e-12, ...
    'Payload-bit-weighted predicted BER calculation failed.');
assert(abs(stats.PredictedTargetMetFrameRatio-0.5) < eps && ...
    abs(stats.FallbackQPSKFrameRatio-0.5) < eps, ...
    'Adaptive frame feasibility/fallback ratios failed.');
assert(~stats.ObservedBER_TargetSatisfied && ~stats.ObservedBER_CI95UpperBelowTarget, ...
    'Measured BER target status was incorrectly reported.');
fprintf('  Test 15 passed: weighted BER prediction and separate target metrics.\n');



% Test 16: frame-cluster bootstrap is reproducible and bounded.
clusterErrors = [0; 2; 1; 4; 0; 3; 1; 2];
clusterBits = 1000*ones(size(clusterErrors));
ciA = cluster_bootstrap_ber_ci(clusterErrors,clusterBits,500,0.95,77);
ciB = cluster_bootstrap_ber_ci(clusterErrors,clusterBits,500,0.95,77);
assert(isequal(ciA,ciB) && ciA(1) <= sum(clusterErrors)/sum(clusterBits) && ...
    ciA(2) >= sum(clusterErrors)/sum(clusterBits) && ciA(1) >= 0 && ciA(2) <= 1, ...
    'Frame-cluster bootstrap interval is not reproducible or valid.');
fprintf('  Test 16 passed: reproducible frame-cluster bootstrap CI.\n');

% Test 17: clustered channel outcomes produce nonzero empirical uncertainty.
clusteredErrors = [0; 0; 0; 0; 20; 20; 20; 20];
clusteredCI = cluster_bootstrap_ber_ci(clusteredErrors,1000*ones(8,1),1000,0.95,91);
assert(clusteredCI(2) > clusteredCI(1), ...
    'Frame-cluster bootstrap failed to capture frame-to-frame variability.');
fprintf('  Test 17 passed: frame-level variability captured by bootstrap.\n');

% Test 18: all-zero bootstrap is explicitly degenerate; Wilson remains a separate metric.
zeroCI = cluster_bootstrap_ber_ci(zeros(6,1),1000*ones(6,1),300,0.95,13);
zeroWilson = bit_error_metrics(0,6000,0,6);
assert(all(zeroCI == 0) && zeroWilson.BER_CI95_Upper > 0, ...
    'Zero-error interval behavior is not clearly separated.');
fprintf('  Test 18 passed: zero-error bootstrap and Wilson interval semantics.\n');

% Test 19: reject malformed frame error counts.
invalidCaught = false;
try
    cluster_bootstrap_ber_ci([0; 4],[100; 3],300,0.95,1);
catch
    invalidCaught = true;
end
assert(invalidCaught,'Cluster bootstrap must reject frame errors exceeding payload bits.');
fprintf('  Test 19 passed: bootstrap input validation.\n');

% Test 20: Alamouti encoder follows the standard orthogonal two-slot pattern.
XinformationTest = randn(N,2,3) + 1i*randn(N,2,3);
XencodedTest = alamouti_encode(XinformationTest);
assert(isequal(size(XencodedTest),[N 2 6]) && ...
    max(abs(XencodedTest(:,1,1)-XinformationTest(:,1,1))) < eps && ...
    max(abs(XencodedTest(:,2,1)-XinformationTest(:,2,1))) < eps && ...
    max(abs(XencodedTest(:,1,2)+conj(XinformationTest(:,2,1)))) < eps && ...
    max(abs(XencodedTest(:,2,2)-conj(XinformationTest(:,1,1)))) < eps, ...
    'Alamouti encoder pattern is inconsistent.');
fprintf('  Test 20 passed: Alamouti orthogonal encoding pattern.\n');

% Test 21: perfect-CSI Alamouti decoding over a frequency-selective 2x2 channel.
cfgSTBC = default_config();
tapsSTBC = generate_mimo_channel(cfgSTBC,'rayleigh');
HstbcTrue = fft(tapsSTBC,N,1);
XinfoSTBC = qam_mod(randi([0 1],N*2*3*2,1),4);
XinfoSTBC = reshape(XinfoSTBC,N,2,3);
XstbcTx = alamouti_encode(XinfoSTBC);
txSTBC = ofdm_modulate(XstbcTx,cp);
rxSTBC = apply_mimo_ofdm_channel(txSTBC,tapsSTBC,0);
Ystbc = ofdm_demodulate(rxSTBC,N,cp);
[XstbcRecovered,stbcPostSinr] = alamouti_decode(Ystbc,HstbcTrue,0);
assert(max(abs(XstbcRecovered(:)-XinfoSTBC(:))) < 1e-8 && ...
    all(isinf(stbcPostSinr(:))), ...
    'Noiseless frequency-selective Alamouti recovery failed.');
fprintf('  Test 21 passed: perfect-CSI Alamouti recovery over multipath MIMO.\n');

% Test 22: noisy Alamouti decoding supports per-subcarrier noise and finite output.
rxSTBCNoisy = apply_mimo_ofdm_channel(txSTBC,tapsSTBC,1e-3/N);
YstbcNoisy = ofdm_demodulate(rxSTBCNoisy,N,cp);
[XstbcNoisy,stbcSinrNoisy] = alamouti_decode(YstbcNoisy,HstbcTrue,1e-3*ones(N,1));
assert(isequal(size(XstbcNoisy),size(XinfoSTBC)) && ...
    all(isfinite(real(XstbcNoisy(:)))) && all(isfinite(imag(XstbcNoisy(:)))) && ...
    all(isfinite(stbcSinrNoisy(:))), ...
    'Noisy Alamouti decoder returned invalid dimensions or non-finite values.');
fprintf('  Test 22 passed: noisy Alamouti combining and per-tone noise support.\n');

% Test 23: reject an odd number of received Alamouti time slots.
oddSlotsCaught = false;
try
    alamouti_decode(Ystbc(:,:,1:3),HstbcTrue,0);
catch
    oddSlotsCaught = true;
end
assert(oddSlotsCaught,'Alamouti decoder must reject an odd number of OFDM slots.');
fprintf('  Test 23 passed: Alamouti odd-slot input validation.\n');

% Test 24: reject invalid Alamouti encoder dimensions.
invalidEncodeCaught = false;
try
    alamouti_encode(randn(N,3,2));
catch
    invalidEncodeCaught = true;
end
assert(invalidEncodeCaught,'Alamouti encoder must require exactly two transmit streams.');
fprintf('  Test 24 passed: Alamouti encoder input validation.\n');

% Test 25: joint selector maximizes net spectral efficiency among feasible candidates.
modOrders = [4 16 64];
predictedModes = [1e-2 8e-4 2e-3; 5e-4 2e-4 5e-4];
efficiencyModes = [2.56 5.12 7.68; 1.28 2.56 3.84];
choice = select_mimo_link_mode(predictedModes,modOrders,1e-3,efficiencyModes);
assert(choice.ModeIndex == 1 && choice.ModulationOrder == 16 && ...
    choice.TargetFeasible && abs(choice.SpectralEfficiency_bpsHz-5.12) < 1e-12, ...
    'Joint mode selector did not maximize feasible spectral efficiency.');
fprintf('  Test 25 passed: highest-efficiency target-feasible joint selection.\n');

% Test 26: when no mode meets target, choose the lowest predicted BER and flag fallback.
predictedFallback = [0.03 0.02 0.01; 0.005 0.006 0.007];
choiceFallback = select_mimo_link_mode(predictedFallback,modOrders,1e-3,efficiencyModes);
assert(choiceFallback.ModeIndex == 2 && choiceFallback.ModulationOrder == 4 && ...
    ~choiceFallback.TargetFeasible && abs(choiceFallback.PredictedBER-0.005) < eps, ...
    'Joint mode selector fallback behavior is incorrect.');
fprintf('  Test 26 passed: minimum-BER fallback with explicit target-infeasible flag.\n');

% Test 27: net spectral efficiency is exactly halved for rate-1 Alamouti vs 2-stream SM.
seSpatial = mimo_mode_spectral_efficiency(cfg,16,'spatial-multiplexing');
seAlamouti = mimo_mode_spectral_efficiency(cfg,16,'alamouti-stbc');
assert(abs(seSpatial/seAlamouti-2) < 1e-12, ...
    'Alamouti versus spatial-multiplexing rate factor is incorrect.');
fprintf('  Test 27 passed: spatial/Alamouti spectral-efficiency accounting.\n');

% Test 28: Alamouti post-combining SINR matches sum(|H|^2)/noise variance.
Hcheck = [1 2; 0.5 1i];
HcheckArray = repmat(reshape(Hcheck,1,2,2),N,1,1);
noiseCheck = 0.25;
[~,sinrCheck] = alamouti_decode(complex(zeros(N,2,2)),HcheckArray,noiseCheck);
expectedSinr = sum(abs(Hcheck(:)).^2)/noiseCheck;
assert(max(abs(sinrCheck(:)-expectedSinr)) < 1e-10, ...
    'Alamouti post-combining SINR does not match the analytical expression.');
fprintf('  Test 28 passed: analytical Alamouti post-combining SINR.\n');

% Test 29: end-to-end adaptive mode runner returns a valid metrics row on a tiny sweep.
cfgSmoke = default_config();
cfgSmoke.adaptiveMimo.snrDb = 12;
cfgSmoke.adaptiveMimo.modulations = [4 16 64];
cfgSmoke.adaptiveMimo.numFrames = 4;
cfgSmoke.adaptiveMimo.targetBER = 1e-3;
cfgSmoke.adaptiveMimo.bootstrapSamples = 200;
cfgSmoke.adaptiveMimo.confidenceLevel = 0.95;
cfgSmoke.adaptiveMimo.bootstrapSeed = 1234;
rng(4321,'twister');
smokeResults = run_adaptive_mimo_mode(cfgSmoke,'awgn');
assert(height(smokeResults) == 1 && isfinite(smokeResults.BER(1)) && ...
    smokeResults.BER(1) >= 0 && smokeResults.BER(1) <= 0.5 && ...
    abs(smokeResults.SpatialMultiplexingFrameRatio(1) + ...
    smokeResults.AlamoutiFrameRatio(1)-1) < 1e-12, ...
    'End-to-end adaptive MIMO mode simulation did not return valid metrics.');
fprintf('  Test 29 passed: end-to-end adaptive mode simulation smoke test.\n');

% Test 30: zero spatial-correlation setting reproduces the independent channel.
cfgCorr0 = default_config();
cfgCorr0.channel.spatialCorrelation.enabled = true;
cfgCorr0.channel.spatialCorrelation.txRho = 0;
cfgCorr0.channel.spatialCorrelation.rxRho = 0;
rng(910,'twister'); tapsIndependent = generate_mimo_channel(cfgCorr0,'rayleigh');
rng(910,'twister'); cfgCorr0.channel.spatialCorrelation.enabled = false;
tapsDisabled = generate_mimo_channel(cfgCorr0,'rayleigh');
assert(max(abs(tapsIndependent(:)-tapsDisabled(:))) < 1e-12, ...
    'Zero spatial correlation must reproduce independent fading.');
fprintf('  Test 30 passed: zero-correlation backward compatibility.\n');

% Test 31: empirical receive-side correlation approaches the configured rho.
cfgRxCorr = default_config();
cfgRxCorr.channel.numTaps = 1; cfgRxCorr.channel.pdp = 1;
cfgRxCorr.channel.spatialCorrelation.enabled = true;
cfgRxCorr.channel.spatialCorrelation.txRho = 0;
cfgRxCorr.channel.spatialCorrelation.rxRho = 0.65;
rng(911,'twister'); sampleRx = complex(zeros(5000,2));
for q = 1:size(sampleRx,1)
    hq = generate_mimo_channel(cfgRxCorr,'rayleigh');
    sampleRx(q,:) = [hq(1,1,1),hq(1,2,1)];
end
rhoRxMeasured = mean(sampleRx(:,1).*conj(sampleRx(:,2))) / ...
    sqrt(mean(abs(sampleRx(:,1)).^2)*mean(abs(sampleRx(:,2)).^2));
assert(abs(real(rhoRxMeasured)-0.65) < 0.05, ...
    'Receive-side empirical spatial correlation is inconsistent with its target.');
fprintf('  Test 31 passed: receive-side spatial correlation calibration.\n');

% Test 32: empirical transmit-side correlation approaches the configured rho.
cfgTxCorr = cfgRxCorr;
cfgTxCorr.channel.spatialCorrelation.txRho = 0.55;
cfgTxCorr.channel.spatialCorrelation.rxRho = 0;
rng(912,'twister'); sampleTx = complex(zeros(5000,2));
for q = 1:size(sampleTx,1)
    hq = generate_mimo_channel(cfgTxCorr,'rayleigh');
    sampleTx(q,:) = [hq(1,1,1),hq(1,1,2)];
end
rhoTxMeasured = mean(sampleTx(:,1).*conj(sampleTx(:,2))) / ...
    sqrt(mean(abs(sampleTx(:,1)).^2)*mean(abs(sampleTx(:,2)).^2));
assert(abs(real(rhoTxMeasured)-0.55) < 0.05, ...
    'Transmit-side empirical spatial correlation is inconsistent with its target.');
fprintf('  Test 32 passed: transmit-side spatial correlation calibration.\n');

% Test 33: unit-diagonal correlation preserves each link's average tap power.
cfgPower = cfgRxCorr;
cfgPower.channel.spatialCorrelation.txRho = 0.8;
cfgPower.channel.spatialCorrelation.rxRho = 0.8;
rng(913,'twister'); powerSamples = zeros(3000,2,2);
for q = 1:size(powerSamples,1)
    hq = generate_mimo_channel(cfgPower,'rayleigh');
    powerSamples(q,:,:) = abs(squeeze(hq(1,:,:))).^2;
end
meanPowers = squeeze(mean(powerSamples,1));
assert(max(abs(meanPowers(:)-1)) < 0.06, ...
    'Spatial correlation changed the per-link average power.');
fprintf('  Test 33 passed: per-link power preserved under correlation.\n');

% Test 34: Rician LOS mean remains unchanged while the diffuse component is correlated.
cfgRicianCorr = default_config();
cfgRicianCorr.channel.spatialCorrelation.enabled = true;
cfgRicianCorr.channel.spatialCorrelation.txRho = 0.7;
cfgRicianCorr.channel.spatialCorrelation.rxRho = 0.6;
[meanRicianExpected,~] = channel_estimation_prior(cfgRicianCorr,'rician');
rng(914,'twister'); meanTaps = complex(zeros(cfgRicianCorr.channel.numTaps,2,2));
numRicianSamples = 2500;
for q = 1:numRicianSamples
    meanTaps = meanTaps + generate_mimo_channel(cfgRicianCorr,'rician');
end
meanTaps = meanTaps/numRicianSamples;
meanFrequency = fft(meanTaps,cfgRicianCorr.system.nFFT,1);
assert(max(abs(meanFrequency(:)-meanRicianExpected(:))) < 0.035, ...
    'Rician LOS mean changed after adding diffuse spatial correlation.');
fprintf('  Test 34 passed: Rician LOS mean preserved under diffuse correlation.\n');

% Test 35: reject invalid spatial-correlation coefficients.
invalidRhoCaught = false;
try
    cfgInvalidRho = default_config();
    cfgInvalidRho.channel.spatialCorrelation.enabled = true;
    cfgInvalidRho.channel.spatialCorrelation.txRho = 1.01;
    generate_mimo_channel(cfgInvalidRho,'rayleigh');
catch
    invalidRhoCaught = true;
end
assert(invalidRhoCaught,'Spatial-correlation generator must reject rho outside [0,0.99].');
fprintf('  Test 35 passed: spatial-correlation coefficient validation.\n');

% Test 36: zero spatial correlation makes joint LMMSE equal per-link LMMSE.
cfgJoint = default_config();
cfgJoint.system.nFFT = 32;
cfgJoint.system.cpLength = 8;
cfgJoint.channel.spatialCorrelation.enabled = true;
cfgJoint.channel.spatialCorrelation.txRho = 0;
cfgJoint.channel.spatialCorrelation.rxRho = 0;
[meanJoint,covJoint] = channel_estimation_prior(cfgJoint,'rayleigh');
tapsJoint = generate_mimo_channel(cfgJoint,'rayleigh');
HjointTrue = fft(tapsJoint,cfgJoint.system.nFFT,1);
noiseJoint = 0.4;
rng(936,'twister');
HjointLS = HjointTrue + sqrt(noiseJoint/2)*(randn(size(HjointTrue))+1i*randn(size(HjointTrue)));
Hindependent = channel_estimation_lmmse(HjointLS,cfgJoint.channel.pdp,noiseJoint, ...
    cfgJoint.system.pilotPower,meanJoint,covJoint);
[~,~,RrxJoint,RtxJoint] = spatial_correlation_factors(cfgJoint);
HspatialZero = channel_estimation_joint_lmmse(HjointLS,cfgJoint.channel.pdp,noiseJoint, ...
    cfgJoint.system.pilotPower,meanJoint,covJoint,RtxJoint,RrxJoint);
assert(max(abs(Hindependent(:)-HspatialZero(:))) < 1e-8, ...
    'Joint LMMSE must match per-link LMMSE when spatial correlation is zero.');
fprintf('  Test 36 passed: zero-correlation joint/per-link LMMSE equivalence.\n');

% Test 37: joint LMMSE reduces MSE for a spatially correlated channel.
cfgJoint.channel.spatialCorrelation.txRho = 0.8;
cfgJoint.channel.spatialCorrelation.rxRho = 0.8;
[meanJoint,covJoint] = channel_estimation_prior(cfgJoint,'rayleigh');
[~,~,RrxJoint,RtxJoint] = spatial_correlation_factors(cfgJoint);
jointModel = channel_estimation_joint_lmmse_setup(cfgJoint.channel.pdp,noiseJoint, ...
    cfgJoint.system.pilotPower,cfgJoint.system.nFFT,covJoint,RtxJoint,RrxJoint);
rng(937,'twister');
mseIndependent = 0; mseSpatial = 0; numJointTrials = 250;
for q = 1:numJointTrials
    hq = generate_mimo_channel(cfgJoint,'rayleigh');
    Hq = fft(hq,cfgJoint.system.nFFT,1);
    HqLS = Hq + sqrt(noiseJoint/2)*(randn(size(Hq))+1i*randn(size(Hq)));
    HqIndependent = channel_estimation_lmmse(HqLS,cfgJoint.channel.pdp,noiseJoint, ...
        cfgJoint.system.pilotPower,meanJoint,covJoint);
    HqSpatial = channel_estimation_joint_lmmse_apply(HqLS,meanJoint,jointModel);
    mseIndependent = mseIndependent + mean(abs(HqIndependent(:)-Hq(:)).^2);
    mseSpatial = mseSpatial + mean(abs(HqSpatial(:)-Hq(:)).^2);
end
mseIndependent = mseIndependent/numJointTrials;
mseSpatial = mseSpatial/numJointTrials;
assert(mseSpatial < 0.92*mseIndependent, ...
    'Joint spatial-frequency LMMSE did not provide the expected correlated-channel MSE gain.');
fprintf('  Test 37 passed: joint LMMSE improves correlated-channel MSE.\n');

% Test 38: mean-aware joint LMMSE recovers a noiseless correlated Rician channel.
[meanRicianJoint,covRicianJoint] = channel_estimation_prior(cfgJoint,'rician');
tapsRicianJoint = generate_mimo_channel(cfgJoint,'rician');
HRicianJoint = fft(tapsRicianJoint,cfgJoint.system.nFFT,1);
[~,~,RrxJoint,RtxJoint] = spatial_correlation_factors(cfgJoint);
HRicianJointEstimate = channel_estimation_joint_lmmse(HRicianJoint,cfgJoint.channel.pdp,0, ...
    cfgJoint.system.pilotPower,meanRicianJoint,covRicianJoint,RtxJoint,RrxJoint);
assert(max(abs(HRicianJointEstimate(:)-HRicianJoint(:))) < 1e-8, ...
    'Noiseless joint LMMSE did not recover a valid Rician channel.');
fprintf('  Test 38 passed: noiseless mean-aware joint LMMSE Rician recovery.\n');

% Test 39: posterior uncertainty is valid and falls as pilot noise is reduced.
modelNoisy = channel_estimation_joint_lmmse_setup(cfgJoint.channel.pdp,0.8, ...
    cfgJoint.system.pilotPower,cfgJoint.system.nFFT,covJoint,RtxJoint,RrxJoint);
modelCleaner = channel_estimation_joint_lmmse_setup(cfgJoint.channel.pdp,0.2, ...
    cfgJoint.system.pilotPower,cfgJoint.system.nFFT,covJoint,RtxJoint,RrxJoint);
assert(numel(modelNoisy.posteriorErrorVarianceByTone)==cfgJoint.system.nFFT && ...
    all(isfinite(modelNoisy.posteriorErrorVarianceByTone)) && ...
    all(modelNoisy.posteriorErrorVarianceByTone>=0) && ...
    mean(modelCleaner.posteriorErrorVarianceByTone) < mean(modelNoisy.posteriorErrorVarianceByTone), ...
    'Joint LMMSE posterior uncertainty is invalid or does not decrease with noise.');
fprintf('  Test 39 passed: joint LMMSE posterior uncertainty behavior.\n');

% Test 40: reject a non-positive-semidefinite spatial covariance matrix.
invalidJointCovarianceCaught = false;
try
    channel_estimation_joint_lmmse_setup(cfgJoint.channel.pdp,0.2, ...
        cfgJoint.system.pilotPower,cfgJoint.system.nFFT,1,[1 2;2 1],eye(2));
catch
    invalidJointCovarianceCaught = true;
end
assert(invalidJointCovarianceCaught, ...
    'Joint LMMSE setup must reject non-PSD spatial covariance matrices.');
fprintf('  Test 40 passed: joint LMMSE spatial covariance validation.\n');

% Test 41: joint estimator is integrated into the adaptive MIMO end-to-end runner.
cfgSmokeJoint = cfgSmoke;
cfgSmokeJoint.channel.spatialCorrelation.enabled = true;
cfgSmokeJoint.channel.spatialCorrelation.txRho = 0.65;
cfgSmokeJoint.channel.spatialCorrelation.rxRho = 0.65;
cfgSmokeJoint.estimation.spatialAwareLMMSE = true;
rng(941,'twister');
jointSmokeResults = run_adaptive_mimo_mode(cfgSmokeJoint,'rayleigh');
assert(height(jointSmokeResults)==1 && isfinite(jointSmokeResults.BER(1)) && ...
    jointSmokeResults.BER(1)>=0 && jointSmokeResults.BER(1)<=0.5, ...
    'Joint spatial-aware estimator end-to-end integration failed.');
fprintf('  Test 41 passed: adaptive MIMO runner with joint spatial LMMSE.\n');

% Test 42: frame-level seed schedule makes the adaptive run reproducible.
cfgPair = cfgSmokeJoint;
cfgPair.adaptiveMimo.frameSeedBase = 424242;
rng(1,'twister'); pairRunA = run_adaptive_mimo_mode(cfgPair,'rayleigh');
rng(9876,'twister'); pairRunB = run_adaptive_mimo_mode(cfgPair,'rayleigh');
assert(isequal(pairRunA.BER,pairRunB.BER) && ...
    isequal(pairRunA.SpectralEfficiency_bpsHz,pairRunB.SpectralEfficiency_bpsHz) && ...
    isequal(pairRunA.SpatialMultiplexingFrameRatio,pairRunB.SpatialMultiplexingFrameRatio), ...
    'Frame-level seed schedule did not reproduce adaptive MIMO results.');
fprintf('  Test 42 passed: deterministic frame-seed schedule for paired estimator studies.\n');

% Test 43: frame-level mean bootstrap is deterministic and returns a valid interval.
bootstrapValues = [0.2;0.4;0.5;0.9;1.1;1.3];
meanCiA = frame_bootstrap_mean_ci(bootstrapValues,500,0.95,943);
meanCiB = frame_bootstrap_mean_ci(bootstrapValues,500,0.95,943);
assert(isequal(meanCiA,meanCiB) && all(isfinite(meanCiA)) && ...
    meanCiA(1)<=meanCiA(2), ...
    'Frame mean bootstrap CI is not reproducible or valid.');
fprintf('  Test 43 passed: reproducible frame-mean bootstrap CI.\n');

% Test 44: joint LMMSE exposes mean posterior variance and conservative upper bound.
assert(isfield(modelNoisy,'posteriorMeanErrorVarianceByTone') && ...
    numel(modelNoisy.posteriorMeanErrorVarianceByTone)==cfgJoint.system.nFFT && ...
    all(isfinite(modelNoisy.posteriorMeanErrorVarianceByTone)) && ...
    all(modelNoisy.posteriorMeanErrorVarianceByTone>=0) && ...
    all(modelNoisy.posteriorMeanErrorVarianceByTone <= ...
        modelNoisy.posteriorErrorVarianceByTone+1e-12), ...
    'Joint posterior mean variance or conservative bound is invalid.');
fprintf('  Test 44 passed: mean and worst-link posterior variance outputs.\n');

% Test 45: direct paired estimator validation reports empirical NMSE/MSE.
cfgEstimatorSmoke = default_config();
cfgEstimatorSmoke.system.nFFT = 16;
cfgEstimatorSmoke.system.cpLength = 8;
cfgEstimatorSmoke.channel.spatialCorrelation.enabled = true;
cfgEstimatorSmoke.channel.spatialCorrelation.txRho = 0.7;
cfgEstimatorSmoke.channel.spatialCorrelation.rxRho = 0.7;
estimatorSmoke = channel_estimator_mse_pair(cfgEstimatorSmoke,'rayleigh',8, ...
    24,200,0.95,944);
assert(height(estimatorSmoke)==1 && isfinite(estimatorSmoke.PerLink_NMSE) && ...
    isfinite(estimatorSmoke.Joint_NMSE) && ...
    isfinite(estimatorSmoke.PerLink_EmpiricalMSE_PerCoefficient) && ...
    isfinite(estimatorSmoke.Joint_EmpiricalMSE_PerCoefficient) && ...
    estimatorSmoke.NumFrames==24, ...
    'Direct paired estimator-MSE validation did not return valid metrics.');
fprintf('  Test 45 passed: direct paired estimator NMSE/MSE smoke test.\n');

% Test 46: direct paired estimator validation is reproducible for a fixed seed.
estimatorSmokeRepeat = channel_estimator_mse_pair(cfgEstimatorSmoke,'rayleigh',8, ...
    24,200,0.95,944);
assert(isequal(estimatorSmoke.PerLink_NMSE,estimatorSmokeRepeat.PerLink_NMSE) && ...
    isequal(estimatorSmoke.Joint_NMSE,estimatorSmokeRepeat.Joint_NMSE) && ...
    isequal(estimatorSmoke.Paired_MSE_Improvement_CI95_Lower, ...
        estimatorSmokeRepeat.Paired_MSE_Improvement_CI95_Lower), ...
    'Direct paired estimator validation is not reproducible.');
fprintf('  Test 46 passed: paired estimator-MSE validation reproducibility.\n');

% Test 47: paired estimator validation rejects invalid frame counts.
invalidEstimatorRunCaught = false;
try
    channel_estimator_mse_pair(cfgEstimatorSmoke,'rayleigh',8,1,200,0.95,945);
catch
    invalidEstimatorRunCaught = true;
end
assert(invalidEstimatorRunCaught, ...
    'Direct paired estimator validation must reject fewer than two frames.');
fprintf('  Test 47 passed: paired estimator validation input checks.\n');

% Test 48: with IID spatial links, paired estimator NMSE values remain equivalent.
cfgEstimatorIID = cfgEstimatorSmoke;
cfgEstimatorIID.channel.spatialCorrelation.txRho = 0;
cfgEstimatorIID.channel.spatialCorrelation.rxRho = 0;
estimatorIID = channel_estimator_mse_pair(cfgEstimatorIID,'rayleigh',8, ...
    24,200,0.95,946);
assert(abs(estimatorIID.PerLink_NMSE-estimatorIID.Joint_NMSE)<1e-10, ...
    'Joint and per-link estimator NMSE must match for IID spatial links.');
fprintf('  Test 48 passed: end-to-end IID estimator equivalence.\n');

% Test 49: full posterior link covariance is Hermitian PSD and matches diagonal summaries.
assert(isfield(modelNoisy,'posteriorLinkCovarianceByTone'), ...
    'Joint estimator must expose full posterior link covariance.');
CpostByTone = modelNoisy.posteriorLinkCovarianceByTone;
assert(isequal(size(CpostByTone),[cfgJoint.system.nRx*cfgJoint.system.nTx, ...
    cfgJoint.system.nRx*cfgJoint.system.nTx,cfgJoint.system.nFFT]), ...
    'Posterior link covariance dimensions are incorrect.');
for k = 1:cfgJoint.system.nFFT
    Ck = CpostByTone(:,:,k);
    assert(norm(Ck-Ck','fro') < 1e-9 && min(real(eig((Ck+Ck')/2))) >= -1e-9, ...
        'Posterior link covariance must be Hermitian positive semidefinite.');
    d = max(real(diag(Ck)),0);
    assert(abs(max(d)-modelNoisy.posteriorErrorVarianceByTone(k)) < 1e-8 && ...
        abs(mean(d)-modelNoisy.posteriorMeanErrorVarianceByTone(k)) < 1e-8, ...
        'Full posterior covariance diagonal disagrees with posterior variance summaries.');
end
fprintf('  Test 49 passed: full posterior link covariance PSD and diagonal consistency.\n');

% Test 50: under IID spatial links, colored-noise builder reduces to the prior scalar model.
cfgIIDNoise = default_config();
cfgIIDNoise.channel.spatialCorrelation.enabled = true;
cfgIIDNoise.channel.spatialCorrelation.txRho = 0;
cfgIIDNoise.channel.spatialCorrelation.rxRho = 0;
[~,covScaleIID] = channel_estimation_prior(cfgIIDNoise,'rayleigh');
[~,~,RrxIID,RtxIID] = spatial_correlation_factors(cfgIIDNoise);
thermalNoiseCheck = 0.13;
modelIIDNoise = channel_estimation_joint_lmmse_setup(cfgIIDNoise.channel.pdp, ...
    thermalNoiseCheck,cfgIIDNoise.system.pilotPower,cfgIIDNoise.system.nFFT, ...
    covScaleIID,RtxIID,RrxIID);
RcoloredIID = channel_error_effective_noise_covariance(modelIIDNoise,thermalNoiseCheck,1);
for k = 1:cfgIIDNoise.system.nFFT
    expectedVariance = thermalNoiseCheck + cfgIIDNoise.system.nTx* ...
        modelIIDNoise.posteriorMeanErrorVarianceByTone(k);
    assert(norm(RcoloredIID(:,:,k)-expectedVariance*eye(cfgIIDNoise.system.nRx),'fro') < 1e-8, ...
        'IID posterior covariance must reduce to isotropic effective noise.');
end
fprintf('  Test 50 passed: IID covariance-aware noise reduces to scalar effective noise.\n');

% Test 51: colored-noise MMSE equalizer is equivalent to explicit noise whitening.
rng(951,'twister');
nTestFFT = 12; nTestRx = 2; nTestTx = 2; nTestSymbols = 3;
Hcolored = complex(randn(nTestFFT,nTestRx,nTestTx),randn(nTestFFT,nTestRx,nTestTx));
Ycolored = complex(randn(nTestFFT,nTestRx,nTestSymbols),randn(nTestFFT,nTestRx,nTestSymbols));
Rcolored = complex(zeros(nTestRx,nTestRx,nTestFFT));
Ywhitened = complex(zeros(size(Ycolored)));
Hwhitened = complex(zeros(size(Hcolored)));
for k = 1:nTestFFT
    A = [1 0;0.45+0.2i 1.3];
    Rcolored(:,:,k) = A*A' + 0.3*eye(nTestRx);
    L = chol(Rcolored(:,:,k),'lower');
    Hwhitened(k,:,:) = reshape(L\reshape(Hcolored(k,:,:),nTestRx,nTestTx),1,nTestRx,nTestTx);
    for s = 1:nTestSymbols
        Ywhitened(k,:,s) = (L\reshape(Ycolored(k,:,s),nTestRx,1)).';
    end
end
[Xcolored,SinrColored] = mimo_equalize(Ycolored,Hcolored,Rcolored,'mmse');
[Xwhite,SinrWhite] = mimo_equalize(Ywhitened,Hwhitened,1,'mmse');
assert(max(abs(Xcolored(:)-Xwhite(:))) < 1e-9 && ...
    max(abs(SinrColored(:)-SinrWhite(:))) < 1e-8, ...
    'Covariance-aware MMSE equalizer disagrees with explicit whitening.');
fprintf('  Test 51 passed: colored-noise MMSE equalizer matches explicit whitening.\n');

% Test 52: colored-noise Alamouti combiner is equivalent to explicit whitening.
rng(952,'twister');
HalamColored = complex(randn(nTestFFT,nTestRx,2),randn(nTestFFT,nTestRx,2));
YalamColored = complex(randn(nTestFFT,nTestRx,4),randn(nTestFFT,nTestRx,4));
RalamColored = complex(zeros(nTestRx,nTestRx,nTestFFT));
HalamWhite = complex(zeros(size(HalamColored)));
YalamWhite = complex(zeros(size(YalamColored)));
for k = 1:nTestFFT
    A = [1 0;0.35-0.25i 1.1];
    RalamColored(:,:,k) = A*A' + 0.2*eye(nTestRx);
    L = chol(RalamColored(:,:,k),'lower');
    HalamWhite(k,:,:) = reshape(L\reshape(HalamColored(k,:,:),nTestRx,2),1,nTestRx,2);
    for s = 1:size(YalamColored,3)
        YalamWhite(k,:,s) = (L\reshape(YalamColored(k,:,s),nTestRx,1)).';
    end
end
[XalamColored,SalamColored] = alamouti_decode(YalamColored,HalamColored,RalamColored);
[XalamWhite,SalamWhite] = alamouti_decode(YalamWhite,HalamWhite,1);
assert(max(abs(XalamColored(:)-XalamWhite(:))) < 1e-9 && ...
    max(abs(SalamColored(:)-SalamWhite(:))) < 1e-8, ...
    'Covariance-aware Alamouti combiner disagrees with explicit whitening.');
fprintf('  Test 52 passed: colored-noise Alamouti combiner matches explicit whitening.\n');

% Test 53: covariance-aware detectors reject non-Hermitian and non-PSD noise matrices.
badCovarianceCaught = false;
Rbad = repmat(eye(2),1,1,nTestFFT);
Rbad(1,2,1) = 0.5i;
try
    mimo_equalize(Ycolored,Hcolored,Rbad,'mmse');
catch
    badCovarianceCaught = true;
end
assert(badCovarianceCaught,'MMSE detector must reject non-Hermitian covariance.');
badCovarianceCaught = false;
Rbad = repmat(eye(2),1,1,nTestFFT);
Rbad(:,:,1) = [1 2;2 1];
try
    alamouti_decode(YalamColored,HalamColored,Rbad);
catch
    badCovarianceCaught = true;
end
assert(badCovarianceCaught,'Alamouti detector must reject non-PSD covariance.');
fprintf('  Test 53 passed: covariance-aware detector validation.\n');

% Test 54: scalar and vector noise inputs remain backward compatible.
noiseVector = linspace(0.05,0.2,nTestFFT).';
[XscalarA,SscalarA] = mimo_equalize(Ycolored,Hcolored,noiseVector,'mmse');
XscalarB = complex(zeros(size(XscalarA)));
SscalarB = zeros(size(SscalarA));
for k = 1:nTestFFT
    [xk,sk] = mimo_equalize(Ycolored(k:k,:,:),Hcolored(k:k,:,:),noiseVector(k),'mmse');
    XscalarB(k,:,:) = xk;
    SscalarB(k,:,:) = sk;
end
assert(max(abs(XscalarA(:)-XscalarB(:))) < 1e-10 && ...
    max(abs(SscalarA(:)-SscalarB(:))) < 1e-9, ...
    'Scalar/vector noise backward compatibility was broken.');
fprintf('  Test 54 passed: scalar/vector detector-noise backward compatibility.\n');

% Test 55: MATLAB suppresses a trailing singleton dimension, but a single-
% symbol Yfd matrix is still a valid [nFFT x nRx x 1] detector input.
Ysingle = Ycolored(:,:,1);
[Xsingle,Ssingle] = mimo_equalize(Ysingle,Hcolored,noiseVector,'mmse');
XsingleExpected = XscalarA(:,:,1);
SsingleExpected = SscalarA(:,:,1);
assert(max(abs(Xsingle(:)-XsingleExpected(:))) < 1e-10 && ...
    max(abs(Ssingle(:)-SsingleExpected(:))) < 1e-9, ...
    'Single-symbol 2-D representation was rejected or changed detector output.');
fprintf('  Test 55 passed: single-symbol detector input with MATLAB trailing-singleton semantics.\n');



% Test 56: zero-Doppler sum-of-sinusoids channel is constant in time.
cfgDopplerTest = default_config();
cfgDopplerTest.channel.doppler.maxDopplerHz = 0;
[HtvZero,~,~] = generate_time_varying_mimo_channel( ...
    cfgDopplerTest,'rayleigh',128,20e3,11056);
HtvReference = repmat(HtvZero(1,:,:,:),size(HtvZero,1),1,1,1);
assert(max(abs(HtvZero(:)-HtvReference(:))) < 1e-12, ...
    'A zero-Doppler channel must remain constant over time.');
fprintf('  Test 56 passed: zero-Doppler channel time invariance.\n');

% Test 57: time-varying channel sequence has the documented dimensions and finite values.
[HtvShape,~,~] = generate_time_varying_mimo_channel( ...
    cfgDopplerTest,'rayleigh',257,100e3,11057);
assert(isequal(size(HtvShape),[257 cfgDopplerTest.channel.numTaps ...
    cfgDopplerTest.system.nRx cfgDopplerTest.system.nTx]) && ...
    all(isfinite(real(HtvShape(:)))) && all(isfinite(imag(HtvShape(:)))), ...
    'Time-varying tap sequence dimensions or values are invalid.');
fprintf('  Test 57 passed: time-varying channel dimensions and finite values.\n');

% Test 58: explicit seed makes the time-varying channel reproducible.
[HtvA,~,~] = generate_time_varying_mimo_channel( ...
    cfgDopplerTest,'rayleigh',256,100e3,11058);
[HtvB,~,~] = generate_time_varying_mimo_channel( ...
    cfgDopplerTest,'rayleigh',256,100e3,11058);
[HtvC,~,~] = generate_time_varying_mimo_channel( ...
    cfgDopplerTest,'rayleigh',256,100e3,11059);
assert(isequal(HtvA,HtvB) && max(abs(HtvA(:)-HtvC(:))) > 1e-8, ...
    'Time-varying channel seed reproducibility failed.');
fprintf('  Test 58 passed: reproducible time-varying channel seed.\n');

% Test 59: finite-sum diffuse ACF approximates the isotropic-scattering Jakes J0 law.
cfgJakes = default_config();
cfgJakes.channel.numTaps = 1;
cfgJakes.channel.pdp = 1;
cfgJakes.channel.doppler.maxDopplerHz = 150;
cfgJakes.channel.doppler.numOscillators = 128;
[~,jakesComponents,~] = generate_time_varying_mimo_channel( ...
    cfgJakes,'rayleigh',30000,2000,11059);
jakesLags = [0 1 2 4 6];
[jakesEstimate,jakesTheory] = doppler_temporal_metrics( ...
    jakesComponents.diffuse,150,2000,jakesLags);
assert(max(abs(real(jakesEstimate)-jakesTheory)) < 0.12, ...
    'Sum-of-sinusoids ACF is too far from the Jakes J0 reference.');
fprintf('  Test 59 passed: diffuse Doppler autocorrelation matches Jakes reference.\n');

% Test 60: time-averaged Rayleigh tap powers preserve the normalized PDP.
cfgPDP = default_config();
cfgPDP.channel.doppler.maxDopplerHz = 90;
[~,pdpComponents,~] = generate_time_varying_mimo_channel( ...
    cfgPDP,'rayleigh',18000,200e3,11060);
pdpMeasured = zeros(cfgPDP.channel.numTaps,1);
for ell = 1:cfgPDP.channel.numTaps
    values = pdpComponents.diffuse(:,ell,:,:);
    pdpMeasured(ell) = mean(abs(values(:)).^2);
end
assert(max(abs(pdpMeasured-cfgPDP.channel.pdp(:))) < 0.10, ...
    'Time-varying channel does not preserve the configured PDP.');
fprintf('  Test 60 passed: average Rayleigh tap powers preserve PDP.\n');

% Test 61: Rician first-tap LOS-to-diffuse ensemble power ratio matches K-factor.
% K is an ensemble-power ratio. A single short finite-sum-of-sinusoids
% realization can have a noticeably different time-average diffuse power.
% Average independent seeds to estimate the ensemble ratio robustly. Use one
% tap so this statistical test remains inexpensive while preserving the model.
cfgRicianDoppler = default_config();
cfgRicianDoppler.channel.numTaps = 1;
cfgRicianDoppler.channel.pdp = 1;
cfgRicianDoppler.channel.doppler.maxDopplerHz = 100;
nKSamples = 10000;
numKRealizations = 8;
losPowerSum = 0;
diffusePowerSum = 0;
for kRealization = 1:numKRealizations
    [~,ricianComponents,~] = generate_time_varying_mimo_channel( ...
        cfgRicianDoppler,'rician',nKSamples,200e3,11060+kRealization);
    losFirstTap = ricianComponents.los(:,1,:,:);
    diffuseFirstTap = ricianComponents.diffuse(:,1,:,:);
    losPowerSum = losPowerSum + sum(abs(losFirstTap(:)).^2);
    diffusePowerSum = diffusePowerSum + sum(abs(diffuseFirstTap(:)).^2);
end
measuredK = losPowerSum/diffusePowerSum;
expectedK = 10^(cfgRicianDoppler.channel.ricianK_dB/10);
expectedLosPower = expectedK/(expectedK+1);
measuredLosPower = losPowerSum/(nKSamples*cfgRicianDoppler.system.nRx* ...
    cfgRicianDoppler.system.nTx*numKRealizations);
assert(abs(measuredLosPower/expectedLosPower-1) < 1e-10, ...
    'Rician deterministic LOS component has incorrect configured power.');
assert(abs(measuredK/expectedK-1) < 0.20, ...
    'Rician first-tap ensemble LOS/diffuse power ratio does not match K-factor.');
fprintf(['  Test 61 passed: Rician Doppler channel preserves first-tap K-factor ' ...
    '(ensemble estimate %.3f, configured %.3f).\n'],measuredK,expectedK);

% Test 62: time-varying apply reduces to the static FIR channel for constant taps.
rng(11062,'twister');
Ntv = cfgDopplerTest.system.nFFT;
Ctv = cfgDopplerTest.system.cpLength;
txTv = complex(randn(Ntv+Ctv,cfgDopplerTest.system.nTx,1), ...
    randn(Ntv+Ctv,cfgDopplerTest.system.nTx,1));
staticTapsTv = generate_mimo_channel(cfgDopplerTest,'rayleigh');
totalTvSamples = Ntv+Ctv;
constantSequence = repmat(reshape(staticTapsTv,1,size(staticTapsTv,1), ...
    size(staticTapsTv,2),size(staticTapsTv,3)),totalTvSamples,1,1,1);
yStaticTv = apply_mimo_ofdm_channel(txTv,staticTapsTv,0);
yDynamicTv = apply_time_varying_mimo_ofdm_channel(txTv,constantSequence,0);
assert(max(abs(yStaticTv(:)-yDynamicTv(:))) < 1e-11, ...
    'Time-varying channel application disagrees with static FIR for constant taps.');
fprintf('  Test 62 passed: constant time-varying FIR matches static channel application.\n');

% Test 63: nonzero Doppler creates measurable within-symbol OFDM ICI.
numProbeSymbols = 48;
XdopplerProbe = complex(zeros(cfgDopplerTest.system.nFFT, ...
    cfgDopplerTest.system.nTx,numProbeSymbols));
activeDopplerTone = 9;
XdopplerProbe(activeDopplerTone,1,:) = 1;
txDopplerProbe = ofdm_modulate(XdopplerProbe,cfgDopplerTest.system.cpLength);
fdZeroCfg = cfgDopplerTest;
fdZeroCfg.channel.doppler.maxDopplerHz = 0;
fdHighCfg = cfgDopplerTest;
fdHighCfg.channel.doppler.maxDopplerHz = 300;
fsProbe = cfgDopplerTest.system.nFFT*15e3;
[seqZero,~,~] = generate_time_varying_mimo_channel(fdZeroCfg,'rayleigh', ...
    size(txDopplerProbe,1)*numProbeSymbols,fsProbe,11063);
[seqHigh,~,~] = generate_time_varying_mimo_channel(fdHighCfg,'rayleigh', ...
    size(txDopplerProbe,1)*numProbeSymbols,fsProbe,11064);
yZero = apply_time_varying_mimo_ofdm_channel(txDopplerProbe,seqZero,0);
yHigh = apply_time_varying_mimo_ofdm_channel(txDopplerProbe,seqHigh,0);
Yzero = ofdm_demodulate(yZero,cfgDopplerTest.system.nFFT,cfgDopplerTest.system.cpLength);
Yhigh = ofdm_demodulate(yHigh,cfgDopplerTest.system.nFFT,cfgDopplerTest.system.cpLength);
zeroTotal = sum(abs(Yzero(:)).^2);
zeroDesired = sum(sum(abs(reshape(Yzero(activeDopplerTone,:,:), ...
    cfgDopplerTest.system.nRx,numProbeSymbols)).^2));
highTotal = sum(abs(Yhigh(:)).^2);
highDesired = sum(sum(abs(reshape(Yhigh(activeDopplerTone,:,:), ...
    cfgDopplerTest.system.nRx,numProbeSymbols)).^2));
iciZero = max(zeroTotal-zeroDesired,0)/max(zeroTotal,eps);
iciHigh = max(highTotal-highDesired,0)/max(highTotal,eps);
assert(iciHigh > iciZero+1e-5, ...
    'Nonzero Doppler should increase OFDM off-subcarrier leakage.');
fprintf('  Test 63 passed: Doppler-induced OFDM inter-carrier interference.\n');

% Test 64: configured transmit spatial correlation is retained by Doppler fading.
cfgSpatialDoppler = default_config();
cfgSpatialDoppler.channel.spatialCorrelation.enabled = true;
cfgSpatialDoppler.channel.spatialCorrelation.txRho = 0.7;
cfgSpatialDoppler.channel.spatialCorrelation.rxRho = 0;
[HspatialDoppler,~,~] = generate_time_varying_mimo_channel( ...
    cfgSpatialDoppler,'rayleigh',12000,200e3,11064);
hTx1 = HspatialDoppler(:,1,1,1);
hTx2 = HspatialDoppler(:,1,1,2);
rhoTxMeasured = sum(hTx1.*conj(hTx2))/sqrt(sum(abs(hTx1).^2)*sum(abs(hTx2).^2));
assert(abs(real(rhoTxMeasured)-0.7) < 0.10, ...
    'Time-varying channel did not preserve transmit spatial correlation.');
fprintf('  Test 64 passed: Doppler process preserves spatial correlation.\n');

% Test 65: Doppler generator rejects an aliased maximum Doppler frequency.
badDopplerCaught = false;
badNyquistCfg = cfgDopplerTest;
badNyquistCfg.channel.doppler.maxDopplerHz = 600;
try
    generate_time_varying_mimo_channel(badNyquistCfg,'rayleigh',100,1000,11065);
catch
    badDopplerCaught = true;
end
assert(badDopplerCaught,'Doppler frequency at/above Nyquist must be rejected.');
fprintf('  Test 65 passed: Doppler Nyquist validation.\n');

% Test 66: time-varying channel application rejects mismatched sequence length.
lengthMismatchCaught = false;
try
    apply_time_varying_mimo_ofdm_channel(txTv,constantSequence(1:end-1,:,:,:),0);
catch
    lengthMismatchCaught = true;
end
assert(lengthMismatchCaught,'Mismatched time-varying channel length must be rejected.');
fprintf('  Test 66 passed: time-varying channel length validation.\n');

% Test 67: AWGN reference time-varying channel is the identity MIMO impulse response.
cfgAwgnDoppler = default_config();
[awgnSequence,~,~] = generate_time_varying_mimo_channel( ...
    cfgAwgnDoppler,'awgn',cfgAwgnDoppler.system.nFFT+cfgAwgnDoppler.system.cpLength, ...
    1e6,11067);
awgnStaticTaps = reshape(awgnSequence(1,:,:,:),cfgAwgnDoppler.channel.numTaps, ...
    cfgAwgnDoppler.system.nRx,cfgAwgnDoppler.system.nTx);
awgnFirstError = awgnStaticTaps(1,:,:)-reshape(eye(2),1,2,2);
awgnTail = awgnStaticTaps(2:end,:,:);
assert(max(abs(awgnFirstError(:))) < 1e-12 && max(abs(awgnTail(:))) < 1e-12, ...
    'AWGN reference channel must be a first-tap identity matrix.');
fprintf('  Test 67 passed: AWGN time-varying reference channel identity.\n');

% Test 68: Doppler temporal metrics reject invalid lag indices.
invalidLagCaught = false;
try
    doppler_temporal_metrics(HtvA,10,1000,[0 size(HtvA,1)]);
catch
    invalidLagCaught = true;
end
assert(invalidLagCaught,'Temporal ACF helper must reject lags outside the sequence.');
fprintf('  Test 68 passed: temporal-autocorrelation lag validation.\n');

% Test 69: invalid Doppler configuration is rejected by the main config validator.
badDopplerConfigCaught = false;
badCfg = default_config();
badCfg.channel.doppler.numOscillators = 3;
try
    validate_config(badCfg);
catch
    badDopplerConfigCaught = true;
end
assert(badDopplerConfigCaught,'Main config validator must reject too few Doppler oscillators.');
fprintf('  Test 69 passed: Doppler configuration validation.\n');


% Test 70: causal hold-last tracker exactly reconstructs a constant channel.
trackingSymbols = 12;
trackingPilotIndices = [1 5 9];
HconstantTracking = complex(zeros(4,2,2,trackingSymbols));
HconstantTracking(:) = 0.7+0.2i;
YconstantPilots = HconstantTracking(:,:,:,trackingPilotIndices);
[HholdTracking,~] = temporal_channel_tracking_apply( ...
    YconstantPilots,trackingPilotIndices,trackingSymbols,'hold-last',[],[]);
assert(max(abs(HholdTracking(:)-HconstantTracking(:))) < 1e-12, ...
    'Hold-last tracker failed for a constant noiseless channel.');
fprintf('  Test 70 passed: causal hold-last tracking of a constant channel.\n');

% Test 71: linear interpolation/extrapolation exactly tracks an affine complex sequence.
trackingTime = 1:trackingSymbols;
affineSequence = 1.2+0.25*trackingTime + 1i*(0.4-0.1*trackingTime);
HaffineTracking = reshape(affineSequence,1,1,1,trackingSymbols);
linearPilotIndices = [1 5 9];
YlinearPilots = HaffineTracking(:,:,:,linearPilotIndices);
[HlinearTracking,~] = temporal_channel_tracking_apply( ...
    YlinearPilots,linearPilotIndices,trackingSymbols,'linear',[],[]);
assert(max(abs(HlinearTracking(:)-HaffineTracking(:))) < 1e-11, ...
    'Linear tracker failed to reconstruct an affine complex channel sequence.');
fprintf('  Test 71 passed: linear interpolation and extrapolation.\n');

% Test 72: the Jakes-LMMSE posterior variance decreases with lower pilot noise.
trackingPilotIndices = 1:4:trackingSymbols;
setupNoisy = temporal_channel_tracker_setup(trackingPilotIndices,trackingSymbols, ...
    0.1,120,1/12000,1.0);
setupClean = temporal_channel_tracker_setup(trackingPilotIndices,trackingSymbols, ...
    0.01,120,1/12000,1.0);
assert(all(setupNoisy.posteriorVariance >= 0) && ...
    mean(setupClean.posteriorVariance) < mean(setupNoisy.posteriorVariance), ...
    'Jakes-LMMSE posterior variance must decrease with pilot noise.');
fprintf('  Test 72 passed: Jakes-LMMSE posterior variance versus pilot SNR.\n');

% Test 73: a zero-Doppler LMMSE predictor uses all past pilots coherently.
staticValue = 0.35-0.6i;
staticPilotValues = reshape(staticValue+[0.05 -0.02 0.01],1,1,1,3);
staticPilotIndices = [1 5 9];
staticSetup = temporal_channel_tracker_setup(staticPilotIndices,trackingSymbols, ...
    0.02,0,1/12000,1.0);
[HstaticPrediction,staticPosterior] = temporal_channel_tracking_apply( ...
    staticPilotValues,staticPilotIndices,trackingSymbols,'jakes-lmmse',[],staticSetup);
assert(all(isfinite(real(HstaticPrediction(:)))) && ...
    all(staticPosterior >= 0) && staticPosterior(end) < staticPosterior(2), ...
    'Zero-Doppler temporal LMMSE prediction or uncertainty behavior failed.');
fprintf('  Test 73 passed: zero-Doppler temporal LMMSE prediction.\n');

% Test 74: a known deterministic Rician LOS prior is preserved exactly.
losMeanTracking = complex(zeros(2,1,1,trackingSymbols));
losMeanTracking(1,1,1,:) = reshape(exp(1i*2*pi*(0:trackingSymbols-1)/17),1,1,1,trackingSymbols);
losMeanTracking(2,1,1,:) = reshape(0.7*exp(-1i*2*pi*(0:trackingSymbols-1)/23),1,1,1,trackingSymbols);
losIndices = [1 4 7 10];
losPilotObs = losMeanTracking(:,:,:,losIndices);
losSetup = temporal_channel_tracker_setup(losIndices,trackingSymbols, ...
    0,80,1/12000,0);
[HlosTracked,losPosterior] = temporal_channel_tracking_apply( ...
    losPilotObs,losIndices,trackingSymbols,'jakes-lmmse',losMeanTracking,losSetup);
assert(max(abs(HlosTracked(:)-losMeanTracking(:))) < 1e-12 && ...
    max(abs(losPosterior)) < 1e-12, ...
    'Temporal tracker did not preserve the deterministic LOS prior.');
fprintf('  Test 74 passed: deterministic Rician LOS prior preservation.\n');

% Test 75: causal Jakes-LMMSE output cannot depend on future pilot observations.
causalIndices = [1 4 7 10];
causalObs = complex(randn(2,1,1,4),randn(2,1,1,4));
causalSetup = temporal_channel_tracker_setup(causalIndices,trackingSymbols, ...
    0.03,150,1/12000,1.0);
[HcausalA,~] = temporal_channel_tracking_apply( ...
    causalObs,causalIndices,trackingSymbols,'jakes-lmmse',[],causalSetup);
causalObsChanged = causalObs;
causalObsChanged(:,:,:,4) = causalObsChanged(:,:,:,4)+100;
[HcausalB,~] = temporal_channel_tracking_apply( ...
    causalObsChanged,causalIndices,trackingSymbols,'jakes-lmmse',[],causalSetup);
causalDifference = HcausalA(:,:,:,1:9)-HcausalB(:,:,:,1:9);
assert(max(abs(causalDifference(:))) < 1e-12, ...
    'Causal Jakes-LMMSE estimate changed when a future pilot was modified.');
fprintf('  Test 75 passed: causal temporal predictor excludes future pilots.\n');

% Test 76: temporal tracker rejects non-increasing pilot indices.
invalidTrackingIndicesCaught = false;
try
    temporal_channel_tracker_setup([1 5 4],trackingSymbols,0.1,100,1/12000,1);
catch
    invalidTrackingIndicesCaught = true;
end
assert(invalidTrackingIndicesCaught,'Unsorted pilot indices must be rejected.');
fprintf('  Test 76 passed: temporal tracker pilot-index validation.\n');

% Test 77: temporal tracker rejects negative pilot-noise variance.
invalidTrackingNoiseCaught = false;
try
    temporal_channel_tracker_setup([1 5 9],trackingSymbols,-0.01,100,1/12000,1);
catch
    invalidTrackingNoiseCaught = true;
end
assert(invalidTrackingNoiseCaught,'Negative pilot-noise variance must be rejected.');
fprintf('  Test 77 passed: temporal tracker noise-variance validation.\n');

% Test 78: temporal tracker posterior variance is finite and bounded by prior variance.
boundedSetup = temporal_channel_tracker_setup([1 5 9],trackingSymbols, ...
    0.05,90,1/12000,0.8);
assert(all(isfinite(boundedSetup.posteriorVariance)) && ...
    all(boundedSetup.posteriorVariance >= -eps) && ...
    all(boundedSetup.posteriorVariance <= 0.8+1e-10), ...
    'Temporal posterior variance must be finite and bounded by prior variance.');
fprintf('  Test 78 passed: temporal posterior variance bounds.\n');

% Test 79: method parsing is case-insensitive for HOLD-LAST.
[HholdUpper,~] = temporal_channel_tracking_apply( ...
    YconstantPilots,trackingPilotIndices,trackingSymbols,'HOLD-LAST',[],[]);
assert(max(abs(HholdUpper(:)-HholdTracking(:))) < 1e-12, ...
    'HOLD-LAST method name must be accepted case-insensitively.');
fprintf('  Test 79 passed: case-insensitive HOLD-LAST method parsing.\n');

% Test 80: method parsing is case-insensitive for LINEAR.
[HlinearUpper,~] = temporal_channel_tracking_apply( ...
    YlinearPilots,linearPilotIndices,trackingSymbols,'LINEAR',[],[]);
assert(max(abs(HlinearUpper(:)-HlinearTracking(:))) < 1e-12, ...
    'LINEAR method name must be accepted case-insensitively.');
fprintf('  Test 80 passed: case-insensitive LINEAR method parsing.\n');

% Test 81: method parsing is case-insensitive for JAKES-LMMSE.
[HjakesUpper,jakesUpperVariance] = temporal_channel_tracking_apply( ...
    staticPilotValues,staticPilotIndices,trackingSymbols,'JAKES-LMMSE',[],staticSetup);
assert(max(abs(HjakesUpper(:)-HstaticPrediction(:))) < 1e-12 && ...
    max(abs(jakesUpperVariance-staticPosterior)) < 1e-12, ...
    'JAKES-LMMSE method name must be accepted case-insensitively.');
fprintf('  Test 81 passed: case-insensitive JAKES-LMMSE method parsing.\n');

% Test 82: orthogonal TDM pilot schedule contains no collisions and reserves data slots.
pilotSchedule82 = mimo_ofdm_pilot_schedule(24,2,4);
assert(isequal(pilotSchedule82.pilotIndicesByTx{1},1:4:21) && ...
    isequal(pilotSchedule82.pilotIndicesByTx{2},2:4:22) && ...
    pilotSchedule82.numPilotSymbols==12 && pilotSchedule82.numDataSymbols==12 && ...
    numel(unique(pilotSchedule82.allPilotIndices))==pilotSchedule82.numPilotSymbols, ...
    'Orthogonal TDM pilot schedule contains an overlap or incorrect indices.');
fprintf('  Test 82 passed: orthogonal MIMO pilot schedule and overhead accounting.\n');

% Test 83: pilot scheduler rejects overlapping per-transmitter pilot slots.
shortSpacingCaught = false;
try
    mimo_ofdm_pilot_schedule(24,2,1);
catch
    shortSpacingCaught = true;
end
assert(shortSpacingCaught,'Pilot spacing shorter than nTx must be rejected.');
fprintf('  Test 83 passed: pilot schedule rejects overlapping pilot occasions.\n');

% Test 84: mapped pilot/data grid activates one Tx at a time in pilot symbols.
[Xmapped84,bitsMapped84,pilotsMapped84,scheduleMapped84] = ...
    build_mimo_ofdm_pilot_data_grid(16,2,12,4,16,1,8201,8202);
for tx = 1:2
    for q = 1:numel(scheduleMapped84.pilotIndicesByTx{tx})
        symbolIndex = scheduleMapped84.pilotIndicesByTx{tx}(q);
        otherTx = 3-tx;
        assert(max(abs(Xmapped84(:,otherTx,symbolIndex)))<1e-12 && ...
            max(abs(Xmapped84(:,tx,symbolIndex)-pilotsMapped84(:,tx,q)))<1e-12, ...
            'Pilot OFDM symbol is not orthogonal across transmit antennas.');
    end
end
dataGridTx1_84 = Xmapped84(:,1,scheduleMapped84.dataIndices);
dataGridTx2_84 = Xmapped84(:,2,scheduleMapped84.dataIndices);
assert(all(size(bitsMapped84)==[16*4 2 scheduleMapped84.numDataSymbols]) && ...
    all(abs(dataGridTx1_84(:))>0) && all(abs(dataGridTx2_84(:))>0), ...
    'Payload bits or data resource elements were not mapped correctly.');
fprintf('  Test 84 passed: full-band known pilots and QAM payload mapping.\n');

% Test 85: pilot extraction recovers synthetic per-link channel observations exactly.
Npilot85 = 16; NrPilot85 = 2; NtPilot85 = 2; NsPilot85 = 12;
[~,~,pilotRef85,sched85] = build_mimo_ofdm_pilot_data_grid( ...
    Npilot85,NtPilot85,NsPilot85,4,4,1,8501,8502);
Hknown85 = complex(randn(Npilot85,NrPilot85,NtPilot85,NsPilot85), ...
    randn(Npilot85,NrPilot85,NtPilot85,NsPilot85));
Ysynthetic85 = complex(zeros(Npilot85,NrPilot85,NsPilot85));
for tx = 1:NtPilot85
    for q = 1:numel(sched85.pilotIndicesByTx{tx})
        idx = sched85.pilotIndicesByTx{tx}(q);
        for rx = 1:NrPilot85
            Ysynthetic85(:,rx,idx) = Hknown85(:,rx,tx,idx).*pilotRef85(:,tx,q);
        end
    end
end
obs85 = extract_mimo_pilot_observations(Ysynthetic85,pilotRef85,sched85);
for tx = 1:NtPilot85
    for q = 1:numel(sched85.pilotIndicesByTx{tx})
        idx = sched85.pilotIndicesByTx{tx}(q);
        expected = Hknown85(:,:,tx,idx);
        pilotError85 = obs85{tx}(:,:,q)-expected;
        assert(max(abs(pilotError85(:)))<1e-11, ...
            'Pilot LS extraction failed for an orthogonal transmit pilot.');
    end
end
fprintf('  Test 85 passed: mapped pilot extraction recovers the isolated MIMO links.\n');

% Test 86: diagonal response of a static FIR channel equals its N-point FFT.
Ndiag86 = 16; cpDiag86 = 4; nsDiag86 = 3; tapsDiag86 = 3;
constantTaps86 = complex(randn(tapsDiag86,2,2),randn(tapsDiag86,2,2));
sequence86 = repmat(reshape(constantTaps86,1,tapsDiag86,2,2), ...
    (Ndiag86+cpDiag86)*nsDiag86,1,1,1);
Hdiag86 = time_varying_ofdm_diagonal_response(sequence86,Ndiag86,cpDiag86,nsDiag86);
expected86 = fft(constantTaps86,Ndiag86,1);
for ss = 1:nsDiag86
    assert(max(abs(reshape(Hdiag86(:,:,:,ss)-expected86,[],1)))<1e-11, ...
        'Static time-varying diagonal response does not match the FIR FFT.');
end
fprintf('  Test 86 passed: static FIR frequency response from sample-varying taps.\n');

% Test 87: diagonal response rejects insufficient CP and inconsistent sample count.
invalidDiagonalCaught = false;
try
    time_varying_ofdm_diagonal_response(sequence86,Ndiag86,1,nsDiag86);
catch
    invalidDiagonalCaught = true;
end
assert(invalidDiagonalCaught,'Diagonal-response helper must reject CP shorter than channel memory.');
fprintf('  Test 87 passed: diagonal OFDM response dimension and CP validation.\n');

% Test 88: end-to-end mapped-pilot time-varying OFDM frame returns finite metrics.
cfgMapped88 = default_config();
cfgMapped88.system.nFFT = 16; cfgMapped88.system.cpLength = 4;
cfgMapped88.system.nTx = 2; cfgMapped88.system.nRx = 2;
cfgMapped88.channel.numTaps = 3; cfgMapped88.channel.pdp = [0.6;0.3;0.1];
cfgMapped88.channel.doppler.numOscillators = 16;
paramsMapped88.nSymbols = 12; paramsMapped88.pilotSpacing = 4;
paramsMapped88.modulationOrder = 4; paramsMapped88.pilotSNRdB = 20;
paramsMapped88.maxDopplerHz = 100; paramsMapped88.txRho = 0.4; paramsMapped88.rxRho = 0.3;
paramsMapped88.dataSeed = 8801; paramsMapped88.pilotSeed = 8802;
paramsMapped88.channelSeed = 8803; paramsMapped88.noiseSeed = 8804;
frameMapped88 = simulate_pilot_mapped_ofdm_frame(cfgMapped88,'rayleigh',paramsMapped88);
assert(numel(frameMapped88.methods)==5 && ...
    all(isfinite([frameMapped88.methods.NMSE])) && ...
    all(isfinite([frameMapped88.methods.BER])) && ...
    isfinite(frameMapped88.ICI_Ratio) && frameMapped88.NumDataSymbols>0, ...
    'End-to-end mapped-pilot OFDM frame returned invalid metrics.');
fprintf('  Test 88 passed: end-to-end mapped-pilot time-varying MIMO-OFDM smoke test.\n');

% Test 89: zero-Doppler mapped-pilot waveform has negligible within-symbol ICI.
paramsMapped89 = paramsMapped88;
paramsMapped89.maxDopplerHz = 0;
frameMapped89 = simulate_pilot_mapped_ofdm_frame(cfgMapped88,'rayleigh',paramsMapped89);
assert(frameMapped89.ICI_Ratio<1e-10, ...
    'A static channel should not create measurable OFDM ICI with a sufficient CP.');
fprintf('  Test 89 passed: zero-Doppler mapped-pilot OFDM ICI reference.\n');

% Test 90: sample-varying Doppler creates measurable ICI in the mapped-pilot frame.
paramsMapped90 = paramsMapped88;
paramsMapped90.maxDopplerHz = 300;
frameMapped90 = simulate_pilot_mapped_ofdm_frame(cfgMapped88,'rayleigh',paramsMapped90);
assert(frameMapped90.ICI_Ratio>frameMapped89.ICI_Ratio+1e-7, ...
    'Increasing Doppler did not increase measured ICI in the actual OFDM waveform.');
fprintf('  Test 90 passed: mapped-pilot waveform captures Doppler-induced ICI.\n');

% Test 91: pilot overhead is reflected in net spectral efficiency.
paramsMapped91 = paramsMapped88;
paramsMapped91.pilotSpacing = 8;
frameMapped91 = simulate_pilot_mapped_ofdm_frame(cfgMapped88,'rayleigh',paramsMapped91);
assert(frameMapped91.PilotOverhead < frameMapped88.PilotOverhead && ...
    frameMapped91.SpectralEfficiency > frameMapped88.SpectralEfficiency, ...
    'Wider pilot spacing should reduce pilot overhead and raise net spectral efficiency.');
fprintf('  Test 91 passed: pilot overhead is included in net spectral efficiency.\n');



% Test 92: for a static multipath channel, the full frequency-coupling
% matrix is block diagonal and exactly matches the per-subcarrier FIR FFT.
Nici92 = 16; cpIci92 = 4; nsIci92 = 3; nTapsIci92 = 3;
staticTaps92 = complex(randn(nTapsIci92,2,2),randn(nTapsIci92,2,2))/sqrt(6);
sequence92 = repmat(reshape(staticTaps92,1,nTapsIci92,2,2), ...
    (Nici92+cpIci92)*nsIci92,1,1,1);
Xici92 = complex(randn(Nici92,2,nsIci92),randn(Nici92,2,nsIci92));
txIci92 = ofdm_modulate(Xici92,cpIci92);
rxIci92 = apply_time_varying_mimo_ofdm_channel(txIci92,sequence92,0);
Yici92 = ofdm_demodulate(rxIci92,Nici92,cpIci92);
Hdiag92 = time_varying_ofdm_diagonal_response(sequence92,Nici92,cpIci92,nsIci92);
for ss = 1:nsIci92
    Hblock92 = build_time_varying_ofdm_ici_matrix(sequence92,Nici92,cpIci92,ss);
    yPred92 = reshape(Hblock92*reshape(Xici92(:,:,ss),[],1),Nici92,2);
    assert(max(abs(yPred92(:)-reshape(Yici92(:,:,ss),[],1)))<1e-10, ...
        'Static full ICI matrix did not reproduce the OFDM waveform.');
    for rr = 1:2
        rows = (rr-1)*Nici92+(1:Nici92);
        for tt = 1:2
            cols = (tt-1)*Nici92+(1:Nici92);
            expectedBlock = diag(reshape(Hdiag92(:,rr,tt,ss),[],1));
            blockError = Hblock92(rows,cols)-expectedBlock;
            assert(max(abs(blockError(:)))<1e-10, ...
                'Static full-channel spatial block differs from per-tone response.');
        end
    end
end
fprintf('  Test 92 passed: static full ICI matrix matches diagonal OFDM channel.\n');

% Test 93: the full time-varying coupling matrix exactly reconstructs the
% noiseless waveform after CP removal and FFT, including off-diagonal ICI.
rng(9301,'twister');
Nici93 = 16; cpIci93 = 4; nsIci93 = 3; nTapsIci93 = 3;
Xici93 = complex(randn(Nici93,2,nsIci93),randn(Nici93,2,nsIci93));
txIci93 = ofdm_modulate(Xici93,cpIci93);
tapsIci93 = 0.15*complex(randn((Nici93+cpIci93)*nsIci93,nTapsIci93,2,2), ...
    randn((Nici93+cpIci93)*nsIci93,nTapsIci93,2,2));
rxIci93 = apply_time_varying_mimo_ofdm_channel(txIci93,tapsIci93,0);
Yici93 = ofdm_demodulate(rxIci93,Nici93,cpIci93);
for ss = 1:nsIci93
    Hblock93 = build_time_varying_ofdm_ici_matrix(tapsIci93,Nici93,cpIci93,ss);
    yPred93 = Hblock93*reshape(Xici93(:,:,ss),[],1);
    yActual93 = reshape(Yici93(:,:,ss),[],1);
    assert(max(abs(yPred93-yActual93))<1e-10, ...
        'Time-varying full ICI matrix failed exact waveform reconstruction.');
end
fprintf('  Test 93 passed: time-varying full ICI matrix reconstructs waveform exactly.\n');

% Test 94: full-matrix MMSE and independent per-tone MMSE agree for a
% block-diagonal (static-channel) OFDM channel.
Nici94 = 16; cpIci94 = 4; nTapsIci94 = 3;
staticTaps94 = complex(randn(nTapsIci94,2,2),randn(nTapsIci94,2,2))/sqrt(8);
sequence94 = repmat(reshape(staticTaps94,1,nTapsIci94,2,2),Nici94+cpIci94,1,1,1);
Xici94 = complex(randn(Nici94,2),randn(Nici94,2));
txIci94 = ofdm_modulate(reshape(Xici94,Nici94,2,1),cpIci94);
rxIci94 = apply_time_varying_mimo_ofdm_channel(txIci94,sequence94,0);
Yici94 = ofdm_demodulate(rxIci94,Nici94,cpIci94);
Hdiag94 = time_varying_ofdm_diagonal_response(sequence94,Nici94,cpIci94,1);
Harray94 = reshape(Hdiag94,Nici94,2,2);
[Xtone94,~] = mimo_equalize(reshape(Yici94,Nici94,2,1),Harray94,0.05,'mmse');
Hblock94 = build_time_varying_ofdm_ici_matrix(sequence94,Nici94,cpIci94,1);
[Xfull94,~] = mimo_ici_equalize(Yici94,Hblock94,0.05,'mmse');
assert(max(abs(Xtone94(:)-Xfull94(:)))<1e-9, ...
    'Full-matrix MMSE must reduce to per-tone MMSE for static channels.');
fprintf('  Test 94 passed: full ICI MMSE equals one-tap MMSE for static channels.\n');

% Test 95: the oracle full-ICI equalizer recovers a noiseless Doppler-shifted
% identity MIMO channel, whereas the full coupling matrix remains finite.
Nici95 = 16; cpIci95 = 4; nsIci95 = 1; fsIci95 = 240e3; fdIci95 = 1800;
Xici95 = complex(randn(Nici95,2),randn(Nici95,2));
txIci95 = ofdm_modulate(reshape(Xici95,Nici95,2,1),cpIci95);
nSamp95 = Nici95+cpIci95;
time95 = (0:nSamp95-1).'/fsIci95;
phase95 = exp(1i*2*pi*fdIci95*time95);
taps95 = complex(zeros(nSamp95,1,2,2));
taps95(:,1,1,1) = phase95;
taps95(:,1,2,2) = phase95;
rxIci95 = apply_time_varying_mimo_ofdm_channel(txIci95,taps95,0);
Yici95 = ofdm_demodulate(rxIci95,Nici95,cpIci95);
Hblock95 = build_time_varying_ofdm_ici_matrix(taps95,Nici95,cpIci95,1);
[Xfull95,~] = mimo_ici_equalize(Yici95,Hblock95,1e-12,'mmse');
assert(max(abs(Xfull95(:)-Xici95(:)))<1e-8 && all(isfinite(Hblock95(:))), ...
    'Oracle full-ICI equalizer failed noiseless Doppler-shifted recovery.');
fprintf('  Test 95 passed: full-ICI oracle recovers noiseless Doppler waveform.\n');

% Test 96: the full ICI equalizer rejects incompatible matrix dimensions.
badIciDimensionsCaught = false;
try
    mimo_ici_equalize(zeros(16,2),zeros(32,31),0.01,'mmse');
catch
    badIciDimensionsCaught = true;
end
assert(badIciDimensionsCaught,'Full ICI equalizer must reject incompatible dimensions.');
fprintf('  Test 96 passed: full ICI equalizer dimension validation.\n');

% Test 97: full ICI equalizer rejects negative noise variance.
badIciNoiseCaught = false;
try
    mimo_ici_equalize(zeros(16,2),eye(32),-0.01,'mmse');
catch
    badIciNoiseCaught = true;
end
assert(badIciNoiseCaught,'Full ICI equalizer must reject negative noise variance.');
fprintf('  Test 97 passed: full ICI equalizer noise validation.\n');

% Test 98: a sample-varying channel produces off-diagonal frequency coupling.
assert(norm(Hblock95-diag(diag(Hblock95)),'fro')>1e-3, ...
    'Doppler-varying channel did not produce off-diagonal ICI coupling.');
fprintf('  Test 98 passed: Doppler channel produces off-diagonal ICI matrix terms.\n');


% Test 99: generalized Alamouti decoder recovers data when channel matrices
% differ between the two codeword slots (noise-free reference).
rng(1599,'twister');
Nvt99 = 24; Nrvt99 = 2;
Hfirst99 = complex(randn(Nvt99,Nrvt99,2),randn(Nvt99,Nrvt99,2));
Hsecond99 = complex(randn(Nvt99,Nrvt99,2),randn(Nvt99,Nrvt99,2));
info99 = complex(randn(Nvt99,2),randn(Nvt99,2));
Ypair99 = complex(zeros(Nvt99,Nrvt99,2));
for kk = 1:Nvt99
    h1 = reshape(Hfirst99(kk,:,:),Nrvt99,2);
    h2 = reshape(Hsecond99(kk,:,:),Nrvt99,2);
    Ypair99(kk,:,1) = (h1*reshape(info99(kk,:),2,1)).';
    s1 = info99(kk,1); s2 = info99(kk,2);
    Ypair99(kk,:,2) = (h2*[-conj(s2);conj(s1)]).';
end
[Xdecoded99,~] = alamouti_decode_timevarying(Ypair99,Hfirst99,Hsecond99,1e-12);
assert(max(abs(Xdecoded99(:)-info99(:)))<1e-8, ...
    'Time-varying Alamouti decoder failed noiseless recovery.');
fprintf('  Test 99 passed: two-slot Alamouti decoding with a changing channel.\n');

% Test 100: for a static channel, the generalized decoder agrees with the
% project's original orthogonal Alamouti combiner.
Hstatic100 = complex(randn(Nvt99,Nrvt99,2),randn(Nvt99,Nrvt99,2));
info100 = complex(randn(Nvt99,2),randn(Nvt99,2));
encoded100 = alamouti_encode(reshape(info100,Nvt99,2,1));
Ystatic100 = complex(zeros(Nvt99,Nrvt99,2));
for slot = 1:2
    for kk = 1:Nvt99
        Hk = reshape(Hstatic100(kk,:,:),Nrvt99,2);
        Ystatic100(kk,:,slot) = (Hk*reshape(encoded100(kk,:,slot),2,1)).';
    end
end
[Xold100,~] = alamouti_decode(Ystatic100,Hstatic100,1e-12);
[Xnew100,~] = alamouti_decode_timevarying(Ystatic100,Hstatic100,Hstatic100,1e-12);
assert(max(abs(Xold100(:)-Xnew100(:)))<1e-7, ...
    'Generalized Alamouti decoder should reduce to the static-channel combiner.');
fprintf('  Test 100 passed: static-channel Alamouti compatibility.\n');

% Test 101: the post-detection SINR and candidate BER predictions respond
% monotonically to a higher receiver-noise variance.
Hpred101 = complex(randn(16,2,2),randn(16,2,2))/sqrt(2);
sePred101 = [4 8;2 4];
predLowNoise101 = predict_pilot_mapped_mcs(Hpred101,0.01,[4 16],1e-2,sePred101);
predHighNoise101 = predict_pilot_mapped_mcs(Hpred101,0.2,[4 16],1e-2,sePred101);
assert(mean(predLowNoise101.postSinrSpatialMultiplexing(:)) > ...
    mean(predHighNoise101.postSinrSpatialMultiplexing(:)) && ...
    all(isfinite(predLowNoise101.predictedBER(:))) && ...
    all(predLowNoise101.predictedBER(:)>=0 & predLowNoise101.predictedBER(:)<=0.5), ...
    'Pilot-mapped MCS predictor failed noise monotonicity or finite BER bounds.');
fprintf('  Test 101 passed: post-SINR-based MCS prediction.\n');

% Test 102: the selector chooses the highest-rate candidate among those
% satisfying the predicted BER constraint.
predMatrix102 = [5e-3 8e-4; 7e-4 2e-4];
seMatrix102 = [8 16;2 4];
choice102 = select_mimo_link_mode(predMatrix102,[4 16],1e-3,seMatrix102);
assert(choice102.TargetFeasible && choice102.ModeIndex==1 && ...
    choice102.ModulationOrder==16 && choice102.SpectralEfficiency_bpsHz==16, ...
    'Adaptive MCS selector did not maximize feasible spectral efficiency.');
fprintf('  Test 102 passed: target-constrained MIMO/QAM mode selection.\n');

% Test 103: V15 end-to-end mapped-pilot frame smoke test with a compact FFT.
cfgV15 = default_config();
cfgV15.system.nFFT = 16; cfgV15.system.cpLength = 4;
cfgV15.channel.numTaps = 3; cfgV15.channel.pdp = [0.6;0.3;0.1];
cfgV15.channel.doppler.numOscillators = 16;
paramsV15.nSymbols = 12; paramsV15.pilotSpacing = 4;
paramsV15.pilotSNRdB = 18; paramsV15.maxDopplerHz = 0;
paramsV15.txRho = 0; paramsV15.rxRho = 0;
paramsV15.dataSeed = 31001; paramsV15.pilotSeed = 31002;
paramsV15.channelSeed = 31003; paramsV15.noiseSeed = 31004;
paramsV15.targetBER = 1e-2; paramsV15.modulations = [4 16];
paramsV15.estimator = 'jakes-lmmse';
frameV15_103 = simulate_pilot_mapped_adaptive_frame(cfgV15,'rayleigh',paramsV15);
assert(any(strcmp(frameV15_103.choice.Mode,{'spatial-multiplexing','alamouti-stbc'})) && ...
    ismember(frameV15_103.choice.ModulationOrder,[4 16]) && ...
    frameV15_103.adaptive.TotalBits>0 && isfinite(frameV15_103.adaptive.BER) && ...
    isfinite(frameV15_103.fixedSM16QAM.BER) && isfinite(frameV15_103.fixedAlamoutiQPSK.BER), ...
    'V15 end-to-end adaptive frame did not return valid results.');
fprintf('  Test 103 passed: pilot-driven adaptive waveform smoke test.\n');

% Test 104: fixed seeds reproduce the same MCS decision and BER.
frameV15_104 = simulate_pilot_mapped_adaptive_frame(cfgV15,'rayleigh',paramsV15);
assert(strcmp(frameV15_103.choice.Mode,frameV15_104.choice.Mode) && ...
    frameV15_103.choice.ModulationOrder==frameV15_104.choice.ModulationOrder && ...
    abs(frameV15_103.adaptive.BER-frameV15_104.adaptive.BER)<eps, ...
    'V15 frame simulation is not reproducible under fixed seeds.');
fprintf('  Test 104 passed: deterministic adaptive-frame reproducibility.\n');

% Test 105: net spectral efficiency includes actual pilot and CP overhead
% and distinguishes two-stream SM from rate-one Alamouti.
smSE105 = frameV15_103.fixedSM16QAM.NetSpectralEfficiency_bpsHz;
alSE105 = frameV15_103.fixedAlamoutiQPSK.NetSpectralEfficiency_bpsHz;
expectedSM105 = 2*4*(frameV15_103.schedule.numDataSymbols/paramsV15.nSymbols)* ...
    cfgV15.system.nFFT/(cfgV15.system.nFFT+cfgV15.system.cpLength);
expectedAL105 = 1*2*(frameV15_103.schedule.numDataSymbols/paramsV15.nSymbols)* ...
    cfgV15.system.nFFT/(cfgV15.system.nFFT+cfgV15.system.cpLength);
assert(abs(smSE105-expectedSM105)<1e-12 && abs(alSE105-expectedAL105)<1e-12 && smSE105>alSE105, ...
    'V15 net spectral efficiency does not reflect mode and pilot/CP overhead.');
fprintf('  Test 105 passed: mode-aware net spectral efficiency.\n');

% Test 106: the frame runner rejects an odd number of payload OFDM slots,
% which cannot be paired into complete Alamouti codewords.
paramsOdd106 = paramsV15; paramsOdd106.nSymbols = 11;
oddFrameCaught106 = false;
try
    simulate_pilot_mapped_adaptive_frame(cfgV15,'rayleigh',paramsOdd106);
catch ME
    oddFrameCaught106 = strcmp(ME.identifier,'simulate_pilot_mapped_adaptive_frame:OddDataSlots');
end
assert(oddFrameCaught106,'V15 must reject odd payload-slot counts for Alamouti pairing.');
fprintf('  Test 106 passed: incomplete Alamouti frame validation.\n');

% Test 107: the time-varying Alamouti detector validates its noise input.
badAlamoutiNoiseCaught107 = false;
try
    alamouti_decode_timevarying(zeros(8,2,2),ones(8,2,2),ones(8,2,2),-0.1);
catch
    badAlamoutiNoiseCaught107 = true;
end
assert(badAlamoutiNoiseCaught107,'Time-varying Alamouti detector must reject negative noise variance.');
fprintf('  Test 107 passed: time-varying Alamouti noise validation.\n');


% Test 108: temporal tracker accepts a per-pilot observation-noise profile.
noiseProfile108 = [0.01 0.03 0.20];
setupProfile108 = temporal_channel_tracker_setup([1 5 9],12,noiseProfile108,120, ...
    1e-4,1.0);
setupScalar108 = temporal_channel_tracker_setup([1 5 9],12,0.01,120,1e-4,1.0);
assert(numel(setupProfile108.pilotNoiseVarianceByPilot)==3 && ...
    setupProfile108.posteriorVariance(12)>=setupScalar108.posteriorVariance(12), ...
    'Per-pilot noise profile was not reflected in temporal posterior uncertainty.');
fprintf('  Test 108 passed: per-pilot temporal-LMMSE observation-noise profile.\n');

% Test 109: posterior channel uncertainty and ICI disturbance both make the
% candidate BER predictor more conservative.
Hblock109 = complex(zeros(16,2,2,4));
for kk=1:16
    Hblock109(kk,1,1,:)=1;
    Hblock109(kk,2,2,:)=1;
end
se109=[4 8;2 4];
base109=predict_pilot_mapped_mcs_block(Hblock109,0.01,[4 16],1e-2,se109, ...
    zeros(2,4),0);
uncertainty109=predict_pilot_mapped_mcs_block(Hblock109,0.01,[4 16],1e-2,se109, ...
    0.25*ones(2,4),0);
ici109=predict_pilot_mapped_mcs_block(Hblock109,0.01,[4 16],1e-2,se109, ...
    zeros(2,4),0.25);
assert(mean(uncertainty109.predictedBER(:))>mean(base109.predictedBER(:)) && ...
    mean(ici109.predictedBER(:))>mean(base109.predictedBER(:)), ...
    'Uncertainty/ICI must reduce predicted effective SINR and raise BER.');
fprintf('  Test 109 passed: posterior-error and ICI-aware MCS prediction.\n');

% Test 110: causal pilot innovations expose excess pilot contamination.
schedule110=mimo_ofdm_pilot_schedule(8,2,4);
obs110=cell(2,1);
for tx=1:2
    obs110{tx}=complex(ones(16,2,2));
    obs110{tx}(:,:,2)=1.3*ones(16,2);
end
prior110=complex(zeros(16,2,1,8));
Hfloor110=ones(16,2,2,2);
[iciFirst110,detailFirst110]=estimate_pilot_ici_variance(obs110,schedule110,1, ...
    prior110,1e-3,1,0,1e-4,1,15e3,Hfloor110,[]);
[iciSecond110,detailSecond110]=estimate_pilot_ici_variance(obs110,schedule110,2, ...
    prior110,1e-3,1,0,1e-4,1,15e3,Hfloor110,iciFirst110);
assert(iciFirst110==0 && iciSecond110>0.1 && ...
    all(detailSecond110.excessVarianceByTx>0), ...
    'Pilot innovation did not detect deliberately added observation contamination.');
fprintf('  Test 110 passed: causal pilot-innovation ICI estimation.\n');

% Test 111: V15.1 makes one decision for each pilot/data burst and returns
% valid selected MIMO mode and modulation for every burst.
frameV15_1_111=simulate_pilot_mapped_adaptive_frame_v15_1(cfgV15,'rayleigh',paramsV15);
assert(numel(frameV15_1_111.blocks)==3 && ...
    all(ismember(frameV15_1_111.adaptive.ModeByBlock, ...
    {'spatial-multiplexing','alamouti-stbc'})) && ...
    all(ismember(frameV15_1_111.adaptive.ModulationOrderByBlock,[4 16])), ...
    'V15.1 did not return a valid decision for every data burst.');
fprintf('  Test 111 passed: per-pilot-burst adaptive MIMO/QAM decisions.\n');

% Test 112: the V15.1 frame-level selected-mode and MCS data fractions each sum to one.
assert(abs(frameV15_1_111.adaptive.SpatialMultiplexingFrameRatio + ...
    frameV15_1_111.adaptive.AlamoutiFrameRatio-1)<1e-12 && ...
    abs(frameV15_1_111.adaptive.QPSKFrameRatio + ...
    frameV15_1_111.adaptive.QAM16FrameRatio-1)<1e-12 && ...
    frameV15_1_111.adaptive.TargetFeasibleBlockRatio>=0 && ...
    frameV15_1_111.adaptive.TargetFeasibleBlockRatio<=1, ...
    'V15.1 block fractions must be normalized over data symbols.');
fprintf('  Test 112 passed: normalized burst mode/MCS/feasibility fractions.\n');

% Test 113: the new burst runner is deterministic under fixed channel/data/pilot seeds.
frameV15_1_113=simulate_pilot_mapped_adaptive_frame_v15_1(cfgV15,'rayleigh',paramsV15);
assert(isequal(frameV15_1_111.adaptive.ModeByBlock,frameV15_1_113.adaptive.ModeByBlock) && ...
    isequal(frameV15_1_111.adaptive.ModulationOrderByBlock, ...
    frameV15_1_113.adaptive.ModulationOrderByBlock) && ...
    abs(frameV15_1_111.adaptive.BER-frameV15_1_113.adaptive.BER)<eps, ...
    'V15.1 burst-adaptive frame is not reproducible.');
fprintf('  Test 113 passed: deterministic V15.1 burst-adaptive waveform.\n');

% Test 114: Alamouti receiver accepts different noise levels in the two slots.
H1_114=complex(randn(12,2,2),randn(12,2,2));
H2_114=complex(randn(12,2,2),randn(12,2,2));
info114=complex(randn(12,2),randn(12,2));
encoded114=alamouti_encode(reshape(info114,12,2,1));
Y114=complex(zeros(12,2,2));
for ss=1:2
    for kk=1:12
        Hk=reshape(H1_114(kk,:,:),2,2);
        if ss==2, Hk=reshape(H2_114(kk,:,:),2,2); end
        Y114(kk,:,ss)=(Hk*reshape(encoded114(kk,:,ss),2,1)).';
    end
end
[X114,sinr114]=alamouti_decode_timevarying(Y114,H1_114,H2_114,[0.01 0.08]);
assert(isequal(size(X114),[12 2]) && isequal(size(sinr114),[12 2]) && ...
    all(isfinite(X114(:))) && all(isfinite(sinr114(:))), ...
    'Two-slot unequal-noise Alamouti detection returned invalid outputs.');
fprintf('  Test 114 passed: unequal per-slot-noise Alamouti MMSE detection.\n');

% Test 115: Alamouti detector rejects invalid two-slot noise vectors.
badNoiseVector115=false;
try
    alamouti_decode_timevarying(zeros(8,2,2),ones(8,2,2),ones(8,2,2),[0.1 -0.1]);
catch ME
    badNoiseVector115=strcmp(ME.identifier,'alamouti_decode_timevarying:InvalidNoise');
end
assert(badNoiseVector115,'Alamouti decoder must reject negative per-slot noise variance.');
fprintf('  Test 115 passed: per-slot Alamouti noise validation.\n');

% Test 116: block MCS predictor rejects an odd number of data symbols.
badBlock116=false;
try
    predict_pilot_mapped_mcs_block(Hblock109(:,:,:,1:3),0.01,[4 16],1e-2,se109, ...
        zeros(2,3),0);
catch ME
    badBlock116=strcmp(ME.identifier,'predict_pilot_mapped_mcs_block:BlockShape');
end
assert(badBlock116,'Burst predictor must reject odd Alamouti pairing length.');
fprintf('  Test 116 passed: burst predictor block-length validation.\n');

% Test 117: the effective-noise output includes thermal, channel-posterior and ICI terms.
probe117=predict_pilot_mapped_mcs_block(Hblock109,0.02,[4 16],1e-2,se109, ...
    0.1*ones(2,4),0.03);
assert(all(abs(probe117.effectiveNoiseVarianceBySymbol-(0.25*ones(1,4)))<1e-12), ...
    'Effective disturbance must include thermal, posterior and ICI contributions.');
fprintf('  Test 117 passed: effective receive-noise accounting.\n');

% Test 118: each block decision is associated with the corresponding pilot occasion.
assert(isequal(frameV15_1_111.blocks(1).dataIndices,3:4) && ...
    isequal(frameV15_1_111.blocks(2).dataIndices,7:8) && ...
    isequal(frameV15_1_111.blocks(3).dataIndices,11:12), ...
    'V15.1 data burst partition does not align with the orthogonal pilot schedule.');
fprintf('  Test 118 passed: pilot-aligned data-burst partitioning.\n');

% Test 119: thermal noise is divided by pilot power, while LS-domain ICI
% variance is already normalized and must not be divided a second time.
lsNoise119 = pilot_ls_noise_profile(0.08,4,[0.01 0.03]);
assert(isequal(size(lsNoise119),[1 2]) && ...
    max(abs(lsNoise119-[0.03 0.05]))<1e-12, ...
    'Pilot LS noise conversion incorrectly scaled the ICI variance.');
fprintf('  Test 119 passed: consistent thermal/ICI units in pilot LS noise.\n');

% Test 120: pilot LS noise conversion rejects invalid ICI variance.
badIci120=false;
try
    pilot_ls_noise_profile(0.01,2,[0.01 -0.02]);
catch ME
    badIci120=strcmp(ME.identifier,'pilot_ls_noise_profile:ICIVariance');
end
assert(badIci120,'Pilot LS noise conversion must reject negative ICI variance.');
fprintf('  Test 120 passed: pilot LS ICI-noise validation.\n');

% Test 121: burst predictor post-SINRs agree with the actual MMSE detectors.
Hblock121=complex(randn(12,2,2,4),randn(12,2,2,4));
noise121=0.07;
pred121=predict_pilot_mapped_mcs_block(Hblock121,noise121,[4 16],1e-2, ...
    [8 16;2 4],zeros(2,4),0);
for ss=1:4
    Hs121=reshape(Hblock121(:,:,:,ss),12,2,2);
    [~,sinrSM121]=mimo_equalize(zeros(12,2,1),Hs121,noise121,'mmse');
    sinrErrorSM121=abs(pred121.postSinrSpatialMultiplexing(:,:,ss)-sinrSM121(:,:,1));
    assert(max(sinrErrorSM121(:))<1e-10, ...
        'Predicted SM post-SINR differs from the receiver MMSE detector.');
end
[~,sinrAL121]=alamouti_decode_timevarying(zeros(12,2,2), ...
    reshape(Hblock121(:,:,:,1),12,2,2),reshape(Hblock121(:,:,:,2),12,2,2), ...
    [noise121 noise121]);
sinrErrorAL121=abs(pred121.postSinrAlamouti(:,:,1)-sinrAL121);
assert(max(sinrErrorAL121(:))<1e-10, ...
    'Predicted Alamouti post-SINR differs from the two-slot receiver.');
fprintf('  Test 121 passed: predictor/receiver post-SINR consistency.\n');

% Test 122: the receiver estimator assembly uses each burst's own causal
% channel prediction only, without substituting a future-burst estimate.
blockTemplate122=struct('HestimateForDecision',[],'posteriorVarianceByTx',[]);
blocks122=repmat(blockTemplate122,2,1);
blocks122(1).HestimateForDecision=ones(3,2,2,2);
blocks122(1).posteriorVarianceByTx=0.1*ones(2,2);
blocks122(2).HestimateForDecision=2*ones(3,2,2,2);
blocks122(2).posteriorVarianceByTx=0.2*ones(2,2);
[Hcausal122,postCausal122,dataMask122]=assemble_burst_channel_estimates_v15_1( ...
    blocks122,{3:4,7:8},3,2,2,9);
errFirst122=abs(Hcausal122(:,:,:,3:4)-1);
errSecond122=abs(Hcausal122(:,:,:,7:8)-2);
assert(max(errFirst122(:))<1e-12 && max(errSecond122(:))<1e-12 && ...
    all(all(postCausal122(:,3:4)==0.1)) && all(all(postCausal122(:,7:8)==0.2)) && ...
    isequal(find(dataMask122(:)),[3;4;7;8]), ...
    'Causal receiver CSI assembly mixed channel estimates between bursts.');
fprintf('  Test 122 passed: receiver CSI remains burst-causal.\n');

% Test 123: ICI innovation tracking slices the matching Tx from a full
% multi-transmit Rician LOS prior before invoking the single-link tracker.
schedule123=mimo_ofdm_pilot_schedule(8,2,4);
obs123={complex(ones(16,2,2)),complex(1.1*ones(16,2,2))};
prior123=complex(zeros(16,2,2,8));
prior123(:,:,1,:)=0.2;
prior123(:,:,2,:)=-0.15;
Hfloor123=ones(16,2,2,4);
[ici123,detail123]=estimate_pilot_ici_variance(obs123,schedule123,2, ...
    prior123,1e-3,1,120,1e-4,1,15e3,Hfloor123,0.001);
assert(isfinite(ici123) && ici123>=0 && ...
    numel(detail123.excessVarianceByTx)==2 && ...
    all(isfinite(detail123.excessVarianceByTx)), ...
    'Full MIMO prior must be sliced per Tx for pilot-innovation tracking.');
fprintf('  Test 123 passed: per-transmitter LOS-prior slicing in ICI estimation.\n');

% Test 124: CRC-16/CCITT-FALSE standard check vector (ASCII "123456789" -> 0x29B1).
bytes124=uint8('123456789'); bits124=zeros(8*numel(bytes124),1);
for q=1:numel(bytes124), bits124((q-1)*8+(1:8))=bitget(bytes124(q),8:-1:1); end
cw124=crc16_ccitt_append(bits124);
crcExpected124=[0 0 1 0 1 0 0 1 1 0 1 1 0 0 0 1].';
assert(isequal(cw124(end-15:end),crcExpected124) && crc16_ccitt_verify(cw124), ...
 'CRC-16/CCITT-FALSE does not match the standard check vector.');
fprintf('  Test 124 passed: CRC-16/CCITT-FALSE standard check vector.\n');

% Test 125: CRC rejects a one-bit-corrupted protected payload.
corrupt125=cw124; corrupt125(5)=1-corrupt125(5);
assert(~crc16_ccitt_verify(corrupt125),'CRC must detect this single-bit corruption.');
fprintf('  Test 125 passed: CRC detects payload corruption.\n');

% Test 126: terminated convolutional encoder has the expected rate and zero tail.
input126=mod((0:256).',2); coded126=conv_encode_rate12(input126);
assert(numel(coded126)==2*(numel(input126)+2) && all(isfinite(coded126)) && ...
 all(coded126==0|coded126==1),'Rate-1/2 terminated encoder dimensions are invalid.');
fprintf('  Test 126 passed: terminated rate-1/2 convolutional encoder.\n');

% Test 127: Viterbi decoding recovers an error-free terminated codeword.
assert(isequal(viterbi_decode_rate12(coded126),[input126;0;0]), ...
 'Viterbi decoder failed noiseless terminated recovery.');
fprintf('  Test 127 passed: noiseless Viterbi decoding.\n');

% Test 128: the hard-decision Viterbi decoder corrects a deterministic channel error.
corrupt128=coded126; corrupt128([3 29 177 333])=1-corrupt128([3 29 177 333]);
assert(isequal(viterbi_decode_rate12(corrupt128),[input126;0;0]), ...
 'Viterbi decoder did not correct the selected deterministic channel errors.');
fprintf('  Test 128 passed: deterministic Viterbi error correction.\n');

% Test 129: uncoded CRC packet fills exact QAM capacity and decodes noiselessly.
packet129=packet_crc_encode_v15_2(1536,'uncoded-crc',129);
decoded129=packet_crc_decode_v15_2(packet129.EncodedBits,packet129);
assert(numel(packet129.EncodedBits)==1536 && packet129.PayloadLength==1520 && ...
 decoded129.CRCValid && decoded129.PayloadErrors==0 && decoded129.CodedBitErrors==0, ...
 'Uncoded CRC packet round trip failed.');
fprintf('  Test 129 passed: uncoded CRC packet round trip and capacity.\n');

% Test 130: convolutional-CRC packet fills exact capacity and decodes noiselessly.
packet130=packet_crc_encode_v15_2(1536,'conv-crc-1/2',130);
decoded130=packet_crc_decode_v15_2(packet130.EncodedBits,packet130);
assert(numel(packet130.EncodedBits)==1536 && packet130.PayloadLength==750 && ...
 decoded130.CRCValid && decoded130.PayloadErrors==0 && decoded130.CodedBitErrors==0, ...
 'Convolutional-CRC packet round trip failed.');
fprintf('  Test 130 passed: convolutional-CRC packet round trip and capacity.\n');

% Test 131: terminated Viterbi decoding repairs a one-bit coded-channel error
% before the CRC check, so the original transport block is accepted.
packet131=packet_crc_encode_v15_2(1536,'conv-crc-1/2',131);
rx131=packet131.EncodedBits; rx131(37)=1-rx131(37);
decoded131=packet_crc_decode_v15_2(rx131,packet131);
assert(decoded131.CRCValid && decoded131.PayloadErrors==0 && ...
 decoded131.CodedBitErrors==1,'FEC did not recover a single coded-bit error.');
fprintf('  Test 131 passed: FEC correction restores CRC-valid packet.\n');

% Test 132: insufficient/odd capacities are rejected by the FEC packet builder.
failed132=false;
try, packet_crc_encode_v15_2(1535,'conv-crc-1/2',132); catch, failed132=true; end
assert(failed132,'Odd rate-1/2 encoded capacity must be rejected.');
fprintf('  Test 132 passed: convolutional packet capacity validation.\n');

% Test 133: end-to-end V15.2 pilot-mapped adaptive frame returns all six
% strategy/coding profiles with valid packet and goodput accounting.
cfg133=default_config();
params133.nSymbols=12; params133.pilotSpacing=4; params133.pilotSNRdB=28;
params133.maxDopplerHz=0; params133.txRho=0; params133.rxRho=0;
params133.dataSeed=13301; params133.pilotSeed=13302;
params133.channelSeed=13303; params133.noiseSeed=13304;
params133.targetBER=1e-2; params133.modulations=[4 16];
params133.estimator='jakes-lmmse';
frame133=simulate_pilot_mapped_adaptive_frame_v15_2(cfg133,'rayleigh',params133);
assert(numel(frame133.candidates)==6,'V15.2 must return three strategies by two coding profiles.');
for q=1:numel(frame133.candidates)
    cnd=frame133.candidates{q};
    assert(cnd.PacketCount==3 && cnd.TotalBits>0 && cnd.TotalCodedBits>0 && ...
       isfinite(cnd.PacketErrorRate) && cnd.PacketErrorRate>=0 && cnd.PacketErrorRate<=1 && ...
       isfinite(cnd.CRCGoodput_bpsHz) && cnd.CRCGoodput_bpsHz>=0 && ...
       cnd.CRCGoodput_bpsHz<=cnd.NetSpectralEfficiency_bpsHz+1e-12, ...
       'V15.2 end-to-end candidate packet/goodput metrics are inconsistent.');
end
fprintf('  Test 133 passed: end-to-end adaptive/fixed CRC and FEC smoke test.\n');

% Test 134: vectorized time-varying channel application matches a direct
% scalar reference for complex MIMO taps, including FIR start-up boundaries.
rng(134,'twister');
tx134=randn(7,2,3)+1i*randn(7,2,3);
tap134=randn(21,3,2,2)+1i*randn(21,3,2,2);
fast134=apply_time_varying_mimo_ofdm_channel(tx134,tap134,0);
xSer134=complex(zeros(21,2));
for tx=1:2, xSer134(:,tx)=reshape(tx134(:,tx,:),21,1); end
refSer134=complex(zeros(21,2));
for n=1:21
 for ell=1:min(3,n)
  for rx=1:2
   h=reshape(tap134(n,ell,rx,:),1,2);
   refSer134(n,rx)=refSer134(n,rx)+h*xSer134(n-ell+1,:).';
  end
 end
end
ref134=complex(zeros(7,2,3));
for rx=1:2, ref134(:,rx,:)=reshape(refSer134(:,rx),7,1,3); end
assert(max(abs(fast134(:)-ref134(:)))<1e-11, ...
 'Vectorized time-varying channel does not match the scalar reference.');
fprintf('  Test 134 passed: vectorized channel equals scalar FIR reference.\n');



% Test 135: soft QAM LLR signs match hard Gray demapping at the exact constellation.
for M135=[4 16 64]
    bits135=randi([0 1],1200*log2(M135),1);
    symbols135=qam_mod(bits135,M135);
    llr135=qam_soft_demod_llr(symbols135,M135,1e-4);
    assert(isequal(double(llr135<0),bits135) && ...
        isequal(double(llr135<0),qam_demod(symbols135,M135)), ...
        sprintf('Soft QAM demapper bit ordering failed for M=%d.',M135));
end
fprintf('  Test 135 passed: max-log soft QAM LLR mapping and bit order.\n');

% Test 136: soft Viterbi has exact noiseless terminated-code recovery.
input136=mod((0:512).',2);
coded136=conv_encode_rate12(input136);
llr136=20*(1-2*double(coded136));
assert(isequal(viterbi_decode_rate12_soft(llr136),[input136;0;0]), ...
    'Soft Viterbi failed noiseless terminated recovery.');
fprintf('  Test 136 passed: noiseless soft-decision Viterbi.\n');

% Test 137: weak wrong bit evidence is discounted by soft Viterbi, while
% hard Viterbi treats those low-confidence flips as fully reliable errors.
input137=mod((0:299).',2); coded137=conv_encode_rate12(input137);
llr137=5*(1-2*double(coded137));
weakWrong137=1:5:numel(llr137);
llr137(weakWrong137)=-0.12*llr137(weakWrong137);
soft137=viterbi_decode_rate12_soft(llr137);
hard137=viterbi_decode_rate12(double(llr137<0));
assert(isequal(soft137,[input137;0;0]) && ...
    ~isequal(hard137,[input137;0;0]), ...
    'Soft Viterbi did not exploit the low reliability of flipped observations.');
fprintf('  Test 137 passed: soft Viterbi uses bit reliability.\n');

% Test 138: both uncoded-CRC and convolutional-CRC packets decode with soft LLRs.
for mode138={'uncoded-crc','conv-crc-1/2'}
    packet138=packet_crc_encode_v15_2(1536,mode138{1},138);
    llr138=20*(1-2*double(packet138.EncodedBits));
    decoded138=packet_soft_decode_v15_2(llr138,packet138);
    assert(decoded138.CRCValid && decoded138.PayloadErrors==0 && ...
        decoded138.CodedBitErrors==0,'Soft packet decoder failed a noiseless packet.');
end
fprintf('  Test 138 passed: soft packet decoder and CRC round trips.\n');

% Test 139: soft demapper and packet decoder reject invalid noise/LLR inputs.
failed139a=false; failed139b=false; failed139c=false;
try, qam_soft_demod_llr(1,4,0); catch, failed139a=true; end
try, qam_soft_demod_llr([1;1i],16,[0.1;NaN]); catch, failed139b=true; end
packet139=packet_crc_encode_v15_2(128,'conv-crc-1/2',139);
try, packet_soft_decode_v15_2(NaN(numel(packet139.EncodedBits),1),packet139); catch, failed139c=true; end
assert(failed139a && failed139b && failed139c,'Soft decoder validation accepted invalid inputs.');
fprintf('  Test 139 passed: soft-demapper and packet-LLR validation.\n');

% Test 140: LLR reliability scales with equalizer output-noise variance.
sym140=qam_mod([0;1;1;0],4);
llrLowNoise140=qam_soft_demod_llr(sym140,4,0.01);
llrHighNoise140=qam_soft_demod_llr(sym140,4,0.1);
assert(all(abs(llrLowNoise140)>=abs(llrHighNoise140)-1e-12), ...
    'LLR magnitude must not decrease as the assumed noise variance decreases.');
fprintf('  Test 140 passed: soft LLR confidence follows output-noise variance.\n');

% Test 141: end-to-end V15.2.1 frame carries soft FEC metadata and valid packet metrics.
cfg141=default_config();
params141.nSymbols=12; params141.pilotSpacing=4; params141.pilotSNRdB=28;
params141.maxDopplerHz=0; params141.txRho=0; params141.rxRho=0;
params141.dataSeed=14101; params141.pilotSeed=14102;
params141.channelSeed=14103; params141.noiseSeed=14104;
params141.targetBER=1e-2; params141.modulations=[4 16]; params141.estimator='jakes-lmmse';
frame141=simulate_pilot_mapped_adaptive_frame_v15_2(cfg141,'rician',params141);
assert(numel(frame141.candidates)==6,'Soft-decision integration must retain six strategy/coding candidates.');
for q=1:numel(frame141.candidates)
    cnd=frame141.candidates{q};
    decoderOK=strcmp(cnd.FECDecoder,'max-log soft-decision Viterbi');
    if strcmp(cnd.CodingMode,'uncoded-crc')
        decoderOK=strcmp(cnd.FECDecoder,'not used (uncoded packet)');
    end
    assert(decoderOK && strcmp(cnd.QAMDemapper,'max-log soft LLR') && ...
        cnd.PacketCount==3 && cnd.TotalCodedBits>0 && ...
        isfinite(cnd.PayloadBER) && cnd.PayloadBER>=0 && cnd.PayloadBER<=1 && ...
        isfinite(cnd.CRCGoodput_bpsHz) && cnd.CRCGoodput_bpsHz>=0, ...
        'V15.2.1 soft-decision end-to-end output is invalid.');
end
fprintf('  Test 141 passed: end-to-end soft-decision V15.2.1 frame.\n');



% Test 142: joint link controller starts with an empty 2x2x2 ACK/NACK state.
state142=joint_link_controller_state_v15_3();
assert(isequal(size(state142.Successes),[2 2 2]) && ...
    isequal(size(state142.Failures),[2 2 2]) && ...
    all(state142.TotalObservations(:)==0) && state142.TotalBursts==0, ...
    'Joint controller initial state dimensions/counts are invalid.');
fprintf('  Test 142 passed: joint controller state initialization.\n');

% Synthetic high-SINR prediction used for selector and feedback unit tests.
pred143.predictedBER=1e-6*ones(2,2);
pred143.modulations=[4 16];
pred143.postSinrSpatialMultiplexing=1e4*ones(16,2,4);
pred143.postSinrAlamouti=1e4*ones(16,2,2);

% Test 143: with well-trained CRC feedback and all feasible actions, the
% selector maximizes expected CRC-goodput and chooses SM-16QAM uncoded.
state143=joint_link_controller_state_v15_3();
state143.Successes(:)=100; state143.TotalObservations(:)=100;
[choice143,diag143]=select_joint_mimo_mcs_fec_v15_3(pred143,4,16,4,0.10,state143, ...
    struct('priorWeight',4,'fecCodingGainDb',2,'highBERNoGainThreshold',0.12));
assert(choice143.ModeIndex==1 && choice143.ModulationOrder==16 && ...
    strcmp(choice143.CodingMode,'uncoded-crc') && diag143.FeasibleActionCount==8 && ...
    choice143.TargetPERFeasible, ...
    'Joint controller did not select the highest-goodput feasible action.');
fprintf('  Test 143 passed: maximum CRC-goodput selection under PER constraint.\n');

% Test 144: at an extremely poor predicted BER, the FEC gain prior is disabled
% rather than promising a coding gain outside the validated operating range.
pred144=pred143; pred144.predictedBER=0.20*ones(2,2);
pred144.postSinrSpatialMultiplexing=0.01*ones(16,2,4);
pred144.postSinrAlamouti=0.01*ones(16,2,2);
[~,diag144]=select_joint_mimo_mcs_fec_v15_3(pred144,4,16,4,0.10, ...
    joint_link_controller_state_v15_3(),struct());
a144=diag144.Actions;
idx144=find([a144.ModeIndex]==1 & [a144.ModulationIndex]==1 & [a144.CodingIndex]==2,1);
assert(a144(idx144).CodingGainDb==0 && ...
    abs(a144(idx144).PredictedPayloadBER-pred144.predictedBER(1,1))<1e-12, ...
    'High-BER guard incorrectly applies an FEC coding gain.');
fprintf('  Test 144 passed: conservative high-BER FEC prior guard.\n');

% Test 145: in the supported region, the FEC abstraction applies a positive
% effective-SNR gain and predicts lower BER than the uncoded action.
pred145=pred143; pred145.predictedBER=0.03*ones(2,2);
pred145.postSinrSpatialMultiplexing=8*ones(16,2,4);
pred145.postSinrAlamouti=8*ones(16,2,2);
[~,diag145]=select_joint_mimo_mcs_fec_v15_3(pred145,4,16,4,0.10, ...
    joint_link_controller_state_v15_3(),struct('fecCodingGainDb',2));
a145=diag145.Actions;
u145=find([a145.ModeIndex]==1 & [a145.ModulationIndex]==1 & [a145.CodingIndex]==1,1);
f145=find([a145.ModeIndex]==1 & [a145.ModulationIndex]==1 & [a145.CodingIndex]==2,1);
assert(a145(f145).CodingGainDb>0 && a145(f145).PredictedPayloadBER < a145(u145).PredictedPayloadBER, ...
    'Supported-region soft-FEC predictor did not reflect coding gain.');
fprintf('  Test 145 passed: supported-region soft-FEC BER abstraction.\n');

% Test 146: one CRC ACK increments only the selected action's success count.
choice146=choice143; state146=joint_link_controller_state_v15_3();
state146=update_joint_link_controller_state_v15_3(state146,choice146,false);
assert(state146.Successes(choice146.ModeIndex,choice146.ModulationIndex,choice146.CodingIndex)==1 && ...
    sum(state146.Successes(:))==1 && sum(state146.Failures(:))==0 && state146.TotalBursts==1, ...
    'CRC ACK was not attributed to the selected joint action.');
fprintf('  Test 146 passed: per-action CRC ACK state update.\n');

% Test 147: an ACK lowers the posterior PER estimate for the same action.
state147=joint_link_controller_state_v15_3();
[choice147a,~]=select_joint_mimo_mcs_fec_v15_3(pred145,4,16,4,0.10,state147,struct());
state147=update_joint_link_controller_state_v15_3(state147,choice147a,false);
[choice147b,~]=select_joint_mimo_mcs_fec_v15_3(pred145,4,16,4,0.10,state147,struct());
assert(choice147b.PosteriorPER<choice147a.PosteriorPER, ...
    'A successful CRC outcome did not reduce posterior PER.');
fprintf('  Test 147 passed: CRC ACK lowers posterior packet-error estimate.\n');

% Test 148: a NACK increases the posterior PER estimate for the selected action.
state148=joint_link_controller_state_v15_3();
[choice148a,~]=select_joint_mimo_mcs_fec_v15_3(pred145,4,16,4,0.10,state148,struct());
state148=update_joint_link_controller_state_v15_3(state148,choice148a,true);
[choice148b,~]=select_joint_mimo_mcs_fec_v15_3(pred145,4,16,4,0.10,state148,struct());
assert(choice148b.PosteriorPER>choice148a.PosteriorPER, ...
    'A failed CRC outcome did not increase posterior PER.');
fprintf('  Test 148 passed: CRC NACK increases posterior packet-error estimate.\n');

% Test 149: end-to-end V15.3 appends the joint candidate and makes one
% mode/modulation/coding decision for every pilot-to-pilot burst.
cfg149=default_config();
params149.nSymbols=12; params149.pilotSpacing=4; params149.pilotSNRdB=24;
params149.maxDopplerHz=0; params149.txRho=0; params149.rxRho=0;
params149.dataSeed=14901; params149.pilotSeed=14902;
params149.channelSeed=14903; params149.noiseSeed=14904;
params149.targetBER=1e-2; params149.targetPER=0.15;
params149.modulations=[4 16]; params149.estimator='jakes-lmmse';
frame149=simulate_pilot_mapped_adaptive_frame_v15_3(cfg149,'rician',params149);
joint149=frame149.candidates{7};
assert(numel(frame149.candidates)==7 && numel(frame149.jointBlocks)==3 && ...
    numel(joint149.ModeByBlock)==3 && numel(joint149.ModulationOrderByBlock)==3 && ...
    numel(joint149.CodingModeByBlock)==3 && joint149.PacketCount==3 && ...
    isfinite(joint149.PayloadBER) && isfinite(joint149.PacketErrorRate) && ...
    joint149.CRCGoodput_bpsHz>=0 && joint149.CRCGoodput_bpsHz<=joint149.NetSpectralEfficiency_bpsHz+1e-12, ...
    'V15.3 end-to-end joint candidate is incomplete or has invalid packet metrics.');
fprintf('  Test 149 passed: end-to-end joint MIMO/modulation/FEC frame.\n');

% Test 150: controller ACK/NACK history persists across adjacent frames when
% explicitly passed back by the caller, and only the new bursts are added.
params150=params149; params150.controllerState=frame149.jointControllerState;
params150.dataSeed=15001; params150.pilotSeed=15002;
params150.channelSeed=15003; params150.noiseSeed=15004;
frame150=simulate_pilot_mapped_adaptive_frame_v15_3(cfg149,'rician',params150);
assert(frame150.jointControllerState.TotalBursts==frame149.jointControllerState.TotalBursts+3 && ...
    sum(frame150.jointControllerState.TotalObservations(:))==frame150.jointControllerState.TotalBursts, ...
    'Joint controller feedback state did not persist across frames.');
fprintf('  Test 150 passed: cross-frame causal CRC feedback persistence.\n');

% Test 151: joint controller frame ratios cover every data symbol exactly once.
assert(abs(joint149.SpatialMultiplexingFrameRatio+joint149.AlamoutiFrameRatio-1)<1e-12 && ...
    abs(joint149.QPSKFrameRatio+joint149.QAM16FrameRatio-1)<1e-12 && ...
    abs(joint149.FECFrameRatio+joint149.UncodedFrameRatio-1)<1e-12, ...
    'Joint mode/modulation/FEC frame ratios are not normalized.');
fprintf('  Test 151 passed: normalized joint action fractions.\n');

% Test 152: every selected joint action uses supported modes, modulations and coding.
assert(all(ismember(joint149.ModeByBlock,{'spatial-multiplexing','alamouti-stbc'})) && ...
    all(ismember(joint149.ModulationOrderByBlock,[4 16])) && ...
    all(ismember(joint149.CodingModeByBlock,{'uncoded-crc','conv-crc-1/2'})), ...
    'Joint controller emitted an unsupported link configuration.');
fprintf('  Test 152 passed: joint controller action-domain validation.\n');

% Test 153: fixed seeds and the same initial controller state reproduce the
% selected block actions and end-to-end BER.
frame153=simulate_pilot_mapped_adaptive_frame_v15_3(cfg149,'rician',params149);
assert(isequal(joint149.ModeByBlock,frame153.candidates{7}.ModeByBlock) && ...
    isequal(joint149.ModulationOrderByBlock,frame153.candidates{7}.ModulationOrderByBlock) && ...
    isequal(joint149.CodingModeByBlock,frame153.candidates{7}.CodingModeByBlock) && ...
    abs(joint149.PayloadBER-frame153.candidates{7}.PayloadBER)<eps, ...
    'V15.3 joint frame is not reproducible for fixed seeds/state.');
fprintf('  Test 153 passed: deterministic joint-controller frame.\n');

% Test 154: invalid packet-error targets are rejected by the selector.
failed154=false;
try, select_joint_mimo_mcs_fec_v15_3(pred143,4,16,4,1.0,state143,struct()); catch, failed154=true; end
assert(failed154,'Joint selector accepted targetPER=1.');
fprintf('  Test 154 passed: joint-selector target-PER validation.\n');

% Test 155: malformed feedback state with negative counts is rejected.
badState155=joint_link_controller_state_v15_3(); badState155.Failures(1)= -1;
failed155=false;
try, select_joint_mimo_mcs_fec_v15_3(pred143,4,16,4,0.1,badState155,struct()); catch, failed155=true; end
assert(failed155,'Joint selector accepted negative feedback counts.');
fprintf('  Test 155 passed: joint feedback-state validation.\n');

% Test 156: after a frame, feedback observations equal the number of decoded
% bursts, one CRC outcome per selected action.
assert(sum(frame149.jointControllerState.TotalObservations(:))==numel(frame149.jointBlocks) && ...
    frame149.jointControllerState.TotalBursts==numel(frame149.jointBlocks), ...
    'Joint controller did not record exactly one CRC outcome per burst.');
fprintf('  Test 156 passed: one CRC feedback observation per burst.\n');

% Test 157: the retained V15.2.1 adaptive/fixed references report normalized
% per-frame mode and modulation weights after the V15.3 joint candidate is added.
base157=frame149.candidates{1};
assert(isfinite(base157.SpatialMultiplexingFrameRatio) && ...
    isfinite(base157.AlamoutiFrameRatio) && ...
    abs(base157.SpatialMultiplexingFrameRatio+base157.AlamoutiFrameRatio-1)<1e-12 && ...
    abs(base157.QPSKFrameRatio+base157.QAM16FrameRatio-1)<1e-12 && ...
    isfinite(base157.PredictedBER) && base157.PredictedBER>=0 && base157.PredictedBER<=0.5, ...
    'V15.3 reference candidate frame weighting/ratios are invalid.');
fprintf('  Test 157 passed: legacy reference block-weighting and normalized ratios.\n');

% Test 158: with a sufficiently poor history for every action, the selector
% explicitly uses its fallback and does not claim PER-target feasibility.
state158=joint_link_controller_state_v15_3();
state158.Failures(:)=100; state158.TotalObservations(:)=100;
[choice158,diag158]=select_joint_mimo_mcs_fec_v15_3(pred143,4,16,4,0.10,state158,struct());
assert(diag158.FeasibleActionCount==0 && ~choice158.TargetPERFeasible && ...
    ~isempty(strfind(choice158.SelectionReason,'fallback')), ...
    'No-feasible-action fallback did not preserve the infeasibility flag.');
fprintf('  Test 158 passed: explicit target-infeasible maximum-goodput fallback.\n');

% Test 159: controller state counters must be internally consistent.
badState159=joint_link_controller_state_v15_3();
badState159.Successes(1)=1; badState159.TotalObservations(1)=0;
failed159=false;
try, select_joint_mimo_mcs_fec_v15_3(pred143,4,16,4,0.10,badState159,struct()); catch, failed159=true; end
assert(failed159,'Joint selector accepted inconsistent ACK/NACK observation counts.');
fprintf('  Test 159 passed: feedback-state counter consistency.\n');

% Test 160: the selector rejects SINR arrays whose time/frequency dimensions
% do not match the burst being scheduled.
badPred160=pred143; badPred160.postSinrSpatialMultiplexing=ones(16,2,3);
failed160=false;
try, select_joint_mimo_mcs_fec_v15_3(badPred160,4,16,4,0.10,state143,struct()); catch, failed160=true; end
assert(failed160,'Joint selector accepted a post-SINR array with the wrong burst length.');
fprintf('  Test 160 passed: selector post-SINR dimensional validation.\n');

fprintf('All 160 tests passed successfully.\n\n');
end
