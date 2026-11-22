function [Xhat, postSinr] = alamouti_decode_timevarying(Ypair, Hfirst, Hsecond, noiseVariance)
%ALAMOUTI_DECODE_TIMEVARYING Two-slot Alamouti detector for a changing channel.
% Ypair   : [nFFT x nRx x 2] received frequency-domain OFDM symbols.
% Hfirst  : [nFFT x nRx x 2] channel during the first data slot.
% Hsecond : [nFFT x nRx x 2] channel during the second data slot.
% noiseVariance may be a scalar (same both slots) or [sigma2_slot1 sigma2_slot2].
%
% The augmented observation [y1;conj(y2)] is processed with the actual
% two-slot effective channel. A diagonal augmented noise covariance allows
% the receiver to include different posterior-channel/ICI disturbance levels
% in the two slots. Output is gain-corrected before hard QAM decisions.

if ~isnumeric(Ypair) || ndims(Ypair) > 3 || size(Ypair,3) ~= 2 || isempty(Ypair) || ...
        any(~isfinite(real(Ypair(:)))) || any(~isfinite(imag(Ypair(:))))
    error('alamouti_decode_timevarying:InvalidReceivedSignal', ...
        'Ypair must be finite [nFFT x nRx x 2] data.');
end
validateChannel(Hfirst,'Hfirst');
validateChannel(Hsecond,'Hsecond');
if ~isequal(size(Hfirst),size(Hsecond)) || ...
        size(Hfirst,1) ~= size(Ypair,1) || size(Hfirst,2) ~= size(Ypair,2) || size(Hfirst,3) ~= 2
    error('alamouti_decode_timevarying:DimensionMismatch', ...
        'Ypair and both channel estimates must have compatible [nFFT x nRx x 2] dimensions.');
end
if ~isnumeric(noiseVariance) || ~isreal(noiseVariance) || isempty(noiseVariance) || ...
        ~(isscalar(noiseVariance) || (isvector(noiseVariance) && numel(noiseVariance)==2)) || ...
        any(~isfinite(noiseVariance(:))) || any(noiseVariance(:)<0)
    error('alamouti_decode_timevarying:InvalidNoise', ...
        'noiseVariance must be a nonnegative scalar or a two-element vector.');
end
if isscalar(noiseVariance)
    noiseBySlot = [double(noiseVariance) double(noiseVariance)];
else
    noiseBySlot = double(noiseVariance(:).');
end

[nFFT,nRx,~] = size(Ypair);
Xhat = complex(zeros(nFFT,2));
postSinr = zeros(nFFT,2);
for k = 1:nFFT
    H1 = reshape(Hfirst(k,:,:),nRx,2);
    H2 = reshape(Hsecond(k,:,:),nRx,2);
    A = [H1(:,1), H1(:,2); conj(H2(:,2)), -conj(H2(:,1))];
    y = [reshape(Ypair(k,:,1),nRx,1); conj(reshape(Ypair(k,:,2),nRx,1))];
    diagNoise = [repmat(max(noiseBySlot(1),1e-12),nRx,1); ...
                 repmat(max(noiseBySlot(2),1e-12),nRx,1)];
    R = diag(diagNoise);
    RinvA = R\A;
    Rinv = R\eye(2*nRx);
    W = (A'*RinvA + eye(2))\(A'*Rinv);
    G = W*A;
    gain = diag(G);
    raw = W*y;
    Xhat(k,:) = (raw./(gain+eps)).';
    for stream = 1:2
        desired = abs(G(stream,stream))^2;
        interference = max(sum(abs(G(stream,:)).^2)-desired,0);
        outputNoise = max(real(W(stream,:)*R*W(stream,:)'),0);
        postSinr(k,stream) = desired/max(interference+outputNoise,eps);
    end
end
end

function validateChannel(H,name)
if ~isnumeric(H) || ndims(H) > 3 || isempty(H) || size(H,3) ~= 2 || ...
        any(~isfinite(real(H(:)))) || any(~isfinite(imag(H(:))))
    error('alamouti_decode_timevarying:InvalidChannel', ...
        '%s must be finite [nFFT x nRx x 2] data.',name);
end
end
