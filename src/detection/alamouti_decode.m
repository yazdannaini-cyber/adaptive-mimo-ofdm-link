function [Xhat, postSinr] = alamouti_decode(Yfd, Hest, noiseVariance)
%ALAMOUTI_DECODE Maximum-ratio Alamouti combiner for two transmit antennas.
% Yfd dimensions: [nFFT x nRx x (2*nPairs)].
% Hest dimensions: [nFFT x nRx x 2], assumed constant over each slot pair.
% noiseVariance may be a scalar, an nFFT-element vector, or a receive-noise
% covariance array [nRx x nRx x nFFT]. The covariance form is whitened before
% Alamouti combining and supports spatially colored estimation error.
% Xhat dimensions: [nFFT x 2 x nPairs].

if ~isnumeric(Yfd) || ndims(Yfd) > 3 || isempty(Yfd) || ...
        any(~isfinite(real(Yfd(:)))) || any(~isfinite(imag(Yfd(:))))
    error('alamouti_decode:InvalidReceivedSignal', 'Yfd must be a finite numeric 3-D array.');
end
if ~isnumeric(Hest) || ndims(Hest) > 3 || ...
        size(Hest,1) ~= size(Yfd,1) || size(Hest,2) ~= size(Yfd,2) || size(Hest,3) ~= 2 || ...
        any(~isfinite(real(Hest(:)))) || any(~isfinite(imag(Hest(:))))
    error('alamouti_decode:InvalidChannel', ...
        'Hest must be finite and have dimensions [nFFT x nRx x 2].');
end

[nFFT,nRx,nSlots] = size(Yfd);
if nRx < 1 || mod(nSlots,2) ~= 0
    error('alamouti_decode:InvalidSlots', 'At least one receive antenna and an even number of slots are required.');
end
useCovariance = isnumeric(noiseVariance) && ndims(noiseVariance) == 3;
if useCovariance
    if ~isequal(size(noiseVariance),[nRx nRx nFFT]) || ...
            any(~isfinite(real(noiseVariance(:)))) || any(~isfinite(imag(noiseVariance(:))))
        error('alamouti_decode:NoiseCovarianceSize', ...
            'Noise covariance must have dimensions [nRx x nRx x nFFT] and finite values.');
    end
    noiseByTone = [];
elseif isscalar(noiseVariance)
    noiseByTone = repmat(double(noiseVariance),nFFT,1);
elseif isnumeric(noiseVariance) && numel(noiseVariance) == nFFT
    noiseByTone = double(noiseVariance(:));
else
    error('alamouti_decode:NoiseSize', ...
        'Noise variance must be scalar, an nFFT-vector, or [nRx x nRx x nFFT].');
end
if ~useCovariance && (~isreal(noiseByTone) || any(~isfinite(noiseByTone)) || any(noiseByTone < 0))
    error('alamouti_decode:InvalidNoise', 'Noise variance must be finite and nonnegative.');
end

nPairs = nSlots/2;
Xhat = complex(zeros(nFFT,2,nPairs));
postSinr = zeros(nFFT,1,nPairs);
for k = 1:nFFT
    H = reshape(Hest(k,:,:),nRx,2);
    sigma2 = NaN;
    if useCovariance
        Rv = reshape(noiseVariance(:,:,k),nRx,nRx);
        if norm(Rv-Rv','fro') > 1e-9*max(1,norm(Rv,'fro'))
            error('alamouti_decode:NonHermitianNoiseCovariance', ...
                'Each receive-noise covariance matrix must be Hermitian.');
        end
        Rv = (Rv+Rv')/2;
        [V,D] = eig(Rv);
        lambda = real(diag(D));
        if min(lambda) < -1e-10*max(1,norm(Rv,2))
            error('alamouti_decode:NonPSDCovariance', ...
                'Receive-noise covariance must be positive semidefinite.');
        end
        lambda = max(lambda,0);
        if max(lambda) <= eps
            Hwhite = H;
            whitening = eye(nRx);
            noNoise = true;
        else
            floorValue = max(max(lambda)*1e-12,realmin);
            whitening = diag(1./sqrt(max(lambda,floorValue)))*V';
            Hwhite = whitening*H;
            noNoise = false;
        end
    else
        sigma2 = noiseByTone(k);
        Hwhite = H;
        whitening = eye(nRx);
        noNoise = sigma2 <= eps;
    end

    h1 = Hwhite(:,1);
    h2 = Hwhite(:,2);
    channelGain = sum(abs(h1).^2 + abs(h2).^2);
    denominator = max(channelGain,eps);
    for pairIndex = 1:nPairs
        firstSlot = 2*pairIndex-1;
        secondSlot = 2*pairIndex;
        y1 = whitening*reshape(Yfd(k,:,firstSlot),nRx,1);
        y2 = whitening*reshape(Yfd(k,:,secondSlot),nRx,1);

        % Orthogonal Alamouti combining after noise whitening.
        z1 = sum(conj(h1).*y1 + h2.*conj(y2));
        z2 = sum(conj(h2).*y1 - h1.*conj(y2));
        Xhat(k,1,pairIndex) = z1/denominator;
        Xhat(k,2,pairIndex) = z2/denominator;
        if noNoise
            postSinr(k,1,pairIndex) = Inf;
        elseif useCovariance
            postSinr(k,1,pairIndex) = channelGain;
        else
            postSinr(k,1,pairIndex) = channelGain/sigma2;
        end
    end
end
end
