function efficiency = mimo_mode_spectral_efficiency(cfg, modulationOrder, txMode)
%MIMO_MODE_SPECTRAL_EFFICIENCY Net payload spectral efficiency for a mode.
% Includes time-multiplexed pilot and cyclic-prefix overhead.
% Spatial multiplexing carries nTx symbols/channel use; Alamouti rate is one.
if ~isfield(cfg,'system') || ~all(isfield(cfg.system,{'nFFT','cpLength','nTx','nDataSymbols'}))
    error('mimo_mode_spectral_efficiency:InvalidConfig','Missing required system configuration.');
end
if ~isscalar(modulationOrder) || ~isfinite(modulationOrder) || modulationOrder < 4 || ...
        abs(log2(modulationOrder)-round(log2(modulationOrder))) > eps || ...
        abs(log2(sqrt(modulationOrder))-round(log2(sqrt(modulationOrder)))) > eps
    error('mimo_mode_spectral_efficiency:InvalidModulation','Invalid square QAM order.');
end
N = cfg.system.nFFT;
cp = cfg.system.cpLength;
Nt = cfg.system.nTx;
Ns = cfg.system.nDataSymbols;
if N <= 0 || cp < 0 || cp >= N || Nt < 1 || Ns <= 0 || Ns ~= floor(Ns) || ...
        any([N cp Nt Ns] ~= floor([N cp Nt Ns]))
    error('mimo_mode_spectral_efficiency:InvalidDimensions','Invalid OFDM dimensions.');
end
switch lower(char(txMode))
    case 'spatial-multiplexing'
        streamFactor = Nt;
    case 'alamouti-stbc'
        if Nt ~= 2 || mod(Ns,2) ~= 0
            error('mimo_mode_spectral_efficiency:AlamoutiTopology', ...
                'Alamouti requires two transmit antennas and an even data-symbol count.');
        end
        streamFactor = 1;
    otherwise
        error('mimo_mode_spectral_efficiency:UnknownMode','Unknown MIMO transmission mode.');
end
efficiency = streamFactor*log2(modulationOrder)*Ns/(Ns+Nt)*N/(N+cp);
end
