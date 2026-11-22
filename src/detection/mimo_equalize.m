function [Xhat, postSinr] = mimo_equalize(Yfd, Hest, noiseVariance, method)
%MIMO_EQUALIZE Per-subcarrier ZF or MMSE MIMO detection.
% noiseVariance may be:
%   - nonnegative scalar;
%   - nFFT-element vector of scalar variances;
%   - [nRx x nRx x nFFT] Hermitian receive-noise covariance matrices.
% The matrix form supports colored noise caused by spatially correlated
% channel-estimation error. Xhat dimensions: [nFFT x nTx x nOFDMsymbols].

if ~isnumeric(Yfd) || ndims(Yfd) > 3 || isempty(Yfd) || ...
        any(~isfinite(real(Yfd(:)))) || any(~isfinite(imag(Yfd(:))))
    error('mimo_equalize:InvalidReceivedSignal','Yfd must be a finite numeric 3-D array.');
end
if ~isnumeric(Hest) || ndims(Hest) > 3 || isempty(Hest) || ...
        any(~isfinite(real(Hest(:)))) || any(~isfinite(imag(Hest(:))))
    error('mimo_equalize:InvalidChannel','Hest must be a finite numeric 3-D array.');
end
[nFFT,nRx,nSymbols] = size(Yfd);
[nFFT_H,nRx_H,nTx] = size(Hest);
if nFFT ~= nFFT_H || nRx ~= nRx_H
    error('mimo_equalize:DimensionMismatch','Yfd and Hest dimensions are inconsistent.');
end

useCovariance = isnumeric(noiseVariance) && ndims(noiseVariance) == 3;
if useCovariance
    if ~isequal(size(noiseVariance),[nRx nRx nFFT]) || ...
            any(~isfinite(real(noiseVariance(:)))) || any(~isfinite(imag(noiseVariance(:))))
        error('mimo_equalize:NoiseCovarianceSize', ...
            'Noise covariance must have dimensions [nRx x nRx x nFFT] and finite values.');
    end
    noiseByTone = [];
elseif isscalar(noiseVariance)
    noiseByTone = repmat(double(noiseVariance),nFFT,1);
elseif isnumeric(noiseVariance) && numel(noiseVariance) == nFFT
    noiseByTone = double(noiseVariance(:));
else
    error('mimo_equalize:NoiseSize', ...
        'noiseVariance must be scalar, an nFFT-vector, or [nRx x nRx x nFFT].');
end
if ~useCovariance && (~isreal(noiseByTone) || any(~isfinite(noiseByTone)) || any(noiseByTone < 0))
    error('mimo_equalize:InvalidNoise','Noise variance values must be finite and nonnegative.');
end
method = lower(char(method));
if ~strcmp(method,'zf') && ~strcmp(method,'mmse')
    error('mimo_equalize:UnknownMethod','Method must be ''zf'' or ''mmse''.');
end

Xhat = complex(zeros(nFFT,nTx,nSymbols));
postSinr = zeros(nFFT,nTx,nSymbols);
for k = 1:nFFT
    H = reshape(Hest(k,:,:),nRx,nTx);
    if useCovariance
        Rv = reshape(noiseVariance(:,:,k),nRx,nRx);
        if norm(Rv-Rv','fro') > 1e-9*max(1,norm(Rv,'fro'))
            error('mimo_equalize:NonHermitianNoiseCovariance', ...
                'Each receive-noise covariance matrix must be Hermitian.');
        end
        Rv = (Rv+Rv')/2;
        lambdaMin = min(real(eig(Rv)));
        if lambdaMin < -1e-10*max(1,norm(Rv,2))
            error('mimo_equalize:NonPSDCovariance', ...
                'Receive-noise covariance must be positive semidefinite.');
        end
        if strcmp(method,'zf') || norm(Rv,'fro') <= eps
            W = pinv(H);
        else
            % Regularize only numerically singular PSD covariance matrices,
            % then use linear solves rather than explicitly forming inv(Rv).
            [V,D] = eig(Rv);
            lambda = max(real(diag(D)),0);
            covarianceFloor = max(max(lambda)*1e-12,realmin);
            if min(lambda) < covarianceFloor
                Rv = V*diag(max(lambda,covarianceFloor))*V';
                Rv = (Rv+Rv')/2;
            end
            RinvH = Rv\H;
            Rinv = Rv\eye(nRx);
            W = (H'*RinvH+eye(nTx))\(H'*Rinv);
        end
    else
        sigma2 = noiseByTone(k);
        Rv = sigma2*eye(nRx);
        if strcmp(method,'zf')
            W = pinv(H);
        else
            W = (H'*H+sigma2*eye(nTx))\H';
        end
    end
    G = W*H;
    gain = diag(G);
    for s = 1:nSymbols
        y = reshape(Yfd(k,:,s),nRx,1);
        raw = W*y;
        corrected = raw./(gain+eps);
        Xhat(k,:,s) = corrected.';
        for stream = 1:nTx
            desired = abs(G(stream,stream))^2;
            totalPower = sum(abs(G(stream,:)).^2);
            interference = max(totalPower-desired,0);
            w = W(stream,:);
            outputNoise = max(real(w*Rv*w'),0);
            postSinr(k,stream,s) = desired/max(interference+outputNoise,eps);
        end
    end
end
end
