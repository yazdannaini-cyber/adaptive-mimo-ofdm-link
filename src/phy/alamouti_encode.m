function Xstbc = alamouti_encode(Xinfo)
%ALAMOUTI_ENCODE Rate-1 Alamouti encoder for two transmit antennas.
% Xinfo dimensions: [nFFT x 2 x nPairs], containing (s1,s2) per tone/pair.
% Xstbc dimensions: [nFFT x 2 x (2*nPairs)]. For every subcarrier k:
%   slot 1: [s1, s2]
%   slot 2: [-conj(s2), conj(s1)]
%
% The conventional orthogonal combiner assumes a constant channel over a pair.
% A receiver with per-slot channel estimates (e.g. V15's generalized MMSE
% decoder) can decode this same mapping when the channel changes between slots.

if ~isnumeric(Xinfo) || ndims(Xinfo) > 3 || size(Xinfo,2) ~= 2 || ...
        isempty(Xinfo) || any(~isfinite(real(Xinfo(:)))) || any(~isfinite(imag(Xinfo(:))))
    error('alamouti_encode:InvalidInput', ...
        'Xinfo must be a finite numeric array with dimensions [nFFT x 2 x nPairs].');
end

[nFFT,~,nPairs] = size(Xinfo);
Xstbc = complex(zeros(nFFT,2,2*nPairs));
for pairIndex = 1:nPairs
    s1 = Xinfo(:,1,pairIndex);
    s2 = Xinfo(:,2,pairIndex);
    firstSlot = 2*pairIndex-1;
    secondSlot = 2*pairIndex;
    Xstbc(:,1,firstSlot) = s1;
    Xstbc(:,2,firstSlot) = s2;
    Xstbc(:,1,secondSlot) = -conj(s2);
    Xstbc(:,2,secondSlot) = conj(s1);
end
end
