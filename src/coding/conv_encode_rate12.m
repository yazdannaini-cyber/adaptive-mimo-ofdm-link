function codedBits = conv_encode_rate12(inputBits)
%CONV_ENCODE_RATE12 Terminated rate-1/2 convolutional encoder, K=3, (7,5) octal.
% Two zero tail bits force the trellis to the all-zero terminal state.
inputBits = double(inputBits(:));
if isempty(inputBits) || any(~isfinite(inputBits)) || ...
        any(inputBits ~= 0 & inputBits ~= 1)
    error('conv_encode_rate12:InvalidBits','inputBits must be a nonempty binary vector.');
end
u = [inputBits;0;0];
codedBits = zeros(2*numel(u),1);
state = 0; % [previous input bit, bit before previous] as a two-bit integer
for k = 1:numel(u)
    b = u(k);
    reg = 4*b + 2*bitget(state,2) + bitget(state,1);
    codedBits(2*k-1) = mod(bitget(reg,3)+bitget(reg,2)+bitget(reg,1),2); % 111
    codedBits(2*k)   = mod(bitget(reg,3)+bitget(reg,1),2);              % 101
    state = 2*b + bitget(state,2);
end
end
