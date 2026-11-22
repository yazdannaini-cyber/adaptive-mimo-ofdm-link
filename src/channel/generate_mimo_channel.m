function taps = generate_mimo_channel(cfg, modelName)
%GENERATE_MIMO_CHANNEL Generate a block-static frequency-selective MIMO channel.
% Output dimensions: [nTaps x nRx x nTx].

nTaps = cfg.channel.numTaps;
nRx = cfg.system.nRx;
nTx = cfg.system.nTx;
pdp = cfg.channel.pdp(:);
if numel(pdp) ~= nTaps
    error('generate_mimo_channel:PDPSize', 'PDP length must equal cfg.channel.numTaps.');
end
if cfg.system.cpLength < nTaps-1
    error('generate_mimo_channel:CPTooShort', 'Cyclic prefix must be at least nTaps-1.');
end

modelName = lower(char(modelName));
taps = complex(zeros(nTaps, nRx, nTx));

switch modelName
    case 'awgn'
        if nRx < nTx
            error('generate_mimo_channel:AWGNTopology', 'AWGN reference channel requires nRx >= nTx.');
        end
        for r = 1:nRx
            for t = 1:nTx
                taps(1,r,t) = double(r == t);
            end
        end
    case {'rayleigh','rician'}
        rayleighTaps = complex(zeros(nTaps, nRx, nTx));
        for r = 1:nRx
            for t = 1:nTx
                rayleighTaps(:,r,t) = sqrt(pdp/2) .* ...
                    (randn(nTaps,1) + 1i*randn(nTaps,1));
            end
        end
        % Apply separable Kronecker correlation to the diffuse component.
        % Unit-diagonal Rtx/Rrx preserve each individual link's average PDP.
        [Lrx,Ltx] = spatial_correlation_factors(cfg);
        for ell = 1:nTaps
            Htap = reshape(rayleighTaps(ell,:,:),nRx,nTx);
            rayleighTaps(ell,:,:) = reshape(Lrx*Htap*Ltx',1,nRx,nTx);
        end
        if strcmp(modelName, 'rayleigh')
            taps = rayleighTaps;
        else
            K = 10^(cfg.channel.ricianK_dB/10);
            losTaps = complex(zeros(nTaps, nRx, nTx));
            for r = 1:nRx
                for t = 1:nTx
                    phase = exp(-1i*2*pi*(r-1)*(t-1)/max(nRx,nTx));
                    losTaps(1,r,t) = phase;
                end
            end
            taps = sqrt(K/(K+1))*losTaps + sqrt(1/(K+1))*rayleighTaps;
        end
    otherwise
        error('generate_mimo_channel:UnknownModel', 'Unknown channel model: %s', modelName);
end
end
