function [estimatedACF, theoreticalACF] = doppler_temporal_metrics(diffuseSequence, maxDopplerHz, sampleRateHz, lagSamples)
%DOPPLER_TEMPORAL_METRICS Estimate diffuse fading ACF and compare with Jakes J0.
% diffuseSequence dimensions: [nSamples x nTaps x nRx x nTx].
% The estimate pools normalized correlations over all taps and spatial links.

if ~isnumeric(diffuseSequence) || isempty(diffuseSequence) || ndims(diffuseSequence) > 4 || ...
        any(~isfinite(real(diffuseSequence(:)))) || any(~isfinite(imag(diffuseSequence(:))))
    error('doppler_temporal_metrics:InvalidSequence','diffuseSequence must be finite numeric data.');
end
if ~isscalar(maxDopplerHz) || ~isreal(maxDopplerHz) || ~isfinite(maxDopplerHz) || maxDopplerHz < 0 || ...
        ~isscalar(sampleRateHz) || ~isreal(sampleRateHz) || ~isfinite(sampleRateHz) || sampleRateHz <= 0
    error('doppler_temporal_metrics:InvalidDoppler','Invalid Doppler or sample rate.');
end
if ~isnumeric(lagSamples) || isempty(lagSamples) || any(~isfinite(lagSamples(:))) || ...
        any(lagSamples(:) < 0) || any(lagSamples(:) ~= floor(lagSamples(:))) || ...
        any(lagSamples(:) >= size(diffuseSequence,1))
    error('doppler_temporal_metrics:InvalidLags','lagSamples must be integer lags within the sequence.');
end

lagSamples = lagSamples(:);
estimatedACF = zeros(numel(lagSamples),1);
theoreticalACF = besselj(0,2*pi*maxDopplerHz*lagSamples/sampleRateHz);
for q = 1:numel(lagSamples)
    lag = lagSamples(q);
    numerator = 0;
    denominator = 0;
    for ell = 1:size(diffuseSequence,2)
        for rx = 1:size(diffuseSequence,3)
            for tx = 1:size(diffuseSequence,4)
                h = diffuseSequence(:,ell,rx,tx);
                if lag == 0
                    left = h;
                    right = h;
                else
                    left = h(1:end-lag);
                    right = h(1+lag:end);
                end
                numerator = numerator + sum(conj(left).*right);
                denominator = denominator + sum(abs(left).^2);
            end
        end
    end
    if denominator > 0
        estimatedACF(q) = numerator/denominator;
    else
        estimatedACF(q) = 0;
    end
end
end
