function bits = qam_demod(symbols, M)
%QAM_DEMOD Hard-decision Gray-coded square-QAM demodulator.
% The output bit ordering is the inverse of qam_mod.

k = log2(M);
L = sqrt(M);
m = log2(L);
if abs(k-round(k)) > eps || abs(m-round(m)) > eps || M < 4
    error('qam_demod:InvalidOrder', 'M must be a square power of four (4, 16, 64, ...).');
end

symbols = symbols(:);
normalizer = sqrt((2/3) * (M-1));
z = symbols * normalizer;
indexI = round((real(z) + (L-1))/2);
indexQ = round((imag(z) + (L-1))/2);
indexI = min(max(indexI, 0), L-1);
indexQ = min(max(indexQ, 0), L-1);
grayI = bitxor(uint32(indexI), bitshift(uint32(indexI), -1));
grayQ = bitxor(uint32(indexQ), bitshift(uint32(indexQ), -1));

B = zeros(numel(symbols), k);
for b = 1:m
    B(:,b) = double(bitget(grayI, m-b+1));
    B(:,m+b) = double(bitget(grayQ, m-b+1));
end
bits = reshape(B.', [], 1);
end
