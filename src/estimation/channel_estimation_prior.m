function [meanResponse, covarianceScale] = channel_estimation_prior(cfg, channelName)
%CHANNEL_ESTIMATION_PRIOR Return the channel mean and diffuse covariance scale.
% The prior matches the channel generators in generate_mimo_channel.m:
%   AWGN    : deterministic identity MIMO channel (zero diffuse covariance)
%   Rayleigh: zero-mean diffuse channel with the configured PDP
%   Rician  : deterministic LOS mean plus diffuse channel with variance 1/(K+1)

N = cfg.system.nFFT;
Nr = cfg.system.nRx;
Nt = cfg.system.nTx;
channelName = lower(char(channelName));
meanResponse = complex(zeros(N,Nr,Nt));
covarianceScale = 1;

switch channelName
    case 'awgn'
        % The AWGN reference channel is explicitly defined as H[k] = I.
        covarianceScale = 0;
        for k = 1:N
            for r = 1:Nr
                for t = 1:Nt
                    meanResponse(k,r,t) = double(r == t);
                end
            end
        end
    case 'rayleigh'
        % Zero-mean diffuse channel. The configured PDP is its covariance.
        covarianceScale = 1;
    case 'rician'
        K = 10^(cfg.channel.ricianK_dB/10);
        losTaps = complex(zeros(cfg.channel.numTaps,Nr,Nt));
        for r = 1:Nr
            for t = 1:Nt
                phase = exp(-1i*2*pi*(r-1)*(t-1)/max(Nr,Nt));
                losTaps(1,r,t) = sqrt(K/(K+1))*phase;
            end
        end
        meanResponse = fft(losTaps,N,1);
        covarianceScale = 1/(K+1);
    otherwise
        error('channel_estimation_prior:UnknownChannel', ...
            'Unknown channel model: %s.', channelName);
end
end
