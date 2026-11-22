function Hls = channel_estimation_ls(pilotYfd, pilotPower)
%CHANNEL_ESTIMATION_LS LS estimate using orthogonal time-multiplexed pilots.
% Pilot OFDM symbol t excites only transmit antenna t on every subcarrier.
% pilotYfd dimensions: [nFFT x nRx x nTx].

[nFFT, nRx, nPilotSymbols] = size(pilotYfd);
if nPilotSymbols < 1 || pilotPower <= 0
    error('channel_estimation_ls:InvalidPilot', 'At least one pilot and positive pilot power are required.');
end
Hls = complex(zeros(nFFT, nRx, nPilotSymbols));
for t = 1:nPilotSymbols
    Hls(:,:,t) = pilotYfd(:,:,t) / sqrt(pilotPower);
end
end
