function model = channel_estimation_joint_lmmse_setup(pdp, noiseVariance, pilotPower, nFFT, covarianceScale, Rtx, Rrx)
%CHANNEL_ESTIMATION_JOINT_LMMSE_SETUP Prepare a joint spatial-frequency LMMSE model.
% The channel array is [nFFT x nRx x nTx]. With MATLAB column-major
% vectorization of reshape(H,nFFT,nRx*nTx), the prior covariance is
%   C_h = kron(kron(Rtx,Rrx), R_freq),
% where R_freq is induced by the PDP. The separable eigen-domain solution
% avoids forming or inverting the full (nFFT*nRx*nTx)-square matrix.
%
% The pilot LS noise variance per channel coefficient is noiseVariance/pilotPower.
% The returned posteriorErrorVarianceByTone is the maximum posterior diagonal
% variance over spatial links at each subcarrier (a conservative common bound).
% posteriorLinkCovarianceByTone retains the full posterior covariance among
% spatial links for each subcarrier, enabling covariance-aware detection.

pdp = pdp(:);
if nargin < 5 || isempty(covarianceScale), covarianceScale = 1; end
if nargin < 6 || isempty(Rtx), Rtx = eye(2); end
if nargin < 7 || isempty(Rrx), Rrx = eye(2); end

if ~isnumeric(pdp) || isempty(pdp) || any(~isfinite(pdp)) || any(pdp < 0) || sum(pdp) <= 0 || ...
        ~isscalar(noiseVariance) || ~isreal(noiseVariance) || ~isfinite(noiseVariance) || noiseVariance < 0 || ...
        ~isscalar(pilotPower) || ~isreal(pilotPower) || ~isfinite(pilotPower) || pilotPower <= 0 || ...
        ~isscalar(nFFT) || ~isreal(nFFT) || ~isfinite(nFFT) || nFFT < 1 || nFFT ~= floor(nFFT) || ...
        ~isscalar(covarianceScale) || ~isreal(covarianceScale) || ~isfinite(covarianceScale) || covarianceScale < 0
    error('channel_estimation_joint_lmmse_setup:InvalidScalarInput', ...
        'Invalid PDP, noise, pilot power, FFT size, or covariance scale.');
end
if ~isnumeric(Rtx) || ~ismatrix(Rtx) || size(Rtx,1) ~= size(Rtx,2) || isempty(Rtx) || ...
        ~isnumeric(Rrx) || ~ismatrix(Rrx) || size(Rrx,1) ~= size(Rrx,2) || isempty(Rrx) || ...
        any(~isfinite(real(Rtx(:)))) || any(~isfinite(imag(Rtx(:)))) || ...
        any(~isfinite(real(Rrx(:)))) || any(~isfinite(imag(Rrx(:))))
    error('channel_estimation_joint_lmmse_setup:InvalidCorrelation', ...
        'Rtx and Rrx must be finite, nonempty square numeric matrices.');
end
if norm(Rtx-Rtx','fro') > 1e-9*max(1,norm(Rtx,'fro')) || ...
        norm(Rrx-Rrx','fro') > 1e-9*max(1,norm(Rrx,'fro'))
    error('channel_estimation_joint_lmmse_setup:NonHermitianCorrelation', ...
        'Rtx and Rrx must be Hermitian.');
end

nTx = size(Rtx,1);
nRx = size(Rrx,1);
Rtx = (Rtx+Rtx')/2;
Rrx = (Rrx+Rrx')/2;
[Utx,Dtx] = eig(Rtx);
[Urx,Drx] = eig(Rrx);
lamTx = real(diag(Dtx));
lamRx = real(diag(Drx));
scaleTx = max(1,max(abs(lamTx)));
scaleRx = max(1,max(abs(lamRx)));
if min(lamTx) < -1e-9*scaleTx || min(lamRx) < -1e-9*scaleRx
    error('channel_estimation_joint_lmmse_setup:NonPSDCorrelation', ...
        'Rtx and Rrx must be positive semidefinite.');
end
lamTx = max(lamTx,0);
lamRx = max(lamRx,0);
Rtx = Utx*diag(lamTx)*Utx';
Rrx = Urx*diag(lamRx)*Urx';

k = (0:nFFT-1).';
Rfreq = complex(zeros(nFFT,nFFT));
pdp = pdp/sum(pdp);
for ell = 0:numel(pdp)-1
    steering = exp(-1i*2*pi*k*ell/nFFT);
    Rfreq = Rfreq + pdp(ell+1)*(steering*steering');
end
Rfreq = covarianceScale*(Rfreq+Rfreq')/2;
Rlink = kron(Rtx,Rrx);

[Ufreq,Dfreq] = eig((Rfreq+Rfreq')/2);
[Ulink,Dlink] = eig((Rlink+Rlink')/2);
lamFreq = max(real(diag(Dfreq)),0);
lamLink = max(real(diag(Dlink)),0);
productEigenvalues = lamFreq*lamLink.';
noiseLS = noiseVariance/pilotPower;

if covarianceScale == 0
    gain = zeros(nFFT,nRx*nTx);
    posteriorErrorVarianceByTone = zeros(nFFT,1);
    posteriorMeanErrorVarianceByTone = zeros(nFFT,1);
    posteriorLinkCovarianceByTone = complex(zeros(nRx*nTx,nRx*nTx,nFFT));
    deterministicPrior = true;
elseif noiseLS <= eps
    tolerance = max(productEigenvalues(:))*1e-12;
    if isempty(tolerance) || tolerance <= 0, tolerance = eps; end
    gain = double(productEigenvalues > tolerance);
    posteriorErrorVarianceByTone = zeros(nFFT,1);
    posteriorMeanErrorVarianceByTone = zeros(nFFT,1);
    posteriorLinkCovarianceByTone = complex(zeros(nRx*nTx,nRx*nTx,nFFT));
    deterministicPrior = false;
else
    gain = productEigenvalues ./ (productEigenvalues + noiseLS);
    posteriorEigenvalues = productEigenvalues .* noiseLS ./ (productEigenvalues + noiseLS);
    posteriorDiagonal = (abs(Ufreq).^2) * posteriorEigenvalues * (abs(Ulink).^2).';
    posteriorDiagonal = max(real(posteriorDiagonal),0);
    posteriorMeanErrorVarianceByTone = mean(posteriorDiagonal,2);
    posteriorErrorVarianceByTone = max(posteriorDiagonal,[],2);

    % Preserve the full posterior covariance across spatial links at each tone.
    % This allows the receiver to account for colored channel-estimation error,
    % rather than reducing the uncertainty to one scalar per subcarrier.
    posteriorLinkCovarianceByTone = complex(zeros(nRx*nTx,nRx*nTx,nFFT));
    freqWeights = abs(Ufreq).^2;
    for tone = 1:nFFT
        posteriorLinkEigenvalues = freqWeights(tone,:) * posteriorEigenvalues;
        posteriorLinkCovariance = Ulink * diag(posteriorLinkEigenvalues) * Ulink';
        posteriorLinkCovarianceByTone(:,:,tone) = ...
            (posteriorLinkCovariance+posteriorLinkCovariance')/2;
    end
    deterministicPrior = false;
end

model.nFFT = nFFT;
model.nRx = nRx;
model.nTx = nTx;
model.Ufreq = Ufreq;
model.Ulink = Ulink;
model.gain = gain;
model.posteriorErrorVarianceByTone = posteriorErrorVarianceByTone;
model.posteriorMeanErrorVarianceByTone = posteriorMeanErrorVarianceByTone;
model.posteriorLinkCovarianceByTone = posteriorLinkCovarianceByTone;
model.deterministicPrior = deterministicPrior;
model.noiseLSVariance = noiseLS;
model.covarianceScale = covarianceScale;
model.txCorrelation = Rtx;
model.rxCorrelation = Rrx;
end
