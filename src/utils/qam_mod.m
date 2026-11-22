function symbols = qam_mod(bits, M)
%QAM_MOD Gray-coded square-QAM modulator (4, 16, or 64-QAM).
% Input bits are grouped per symbol. Average constellation energy is 1.

validateattributes(M, {'numeric'}, {'scalar','integer','positive'});
k = log2(M);
L = sqrt(M);
m = log2(L);
if abs(k-round(k)) > eps || abs(m-round(m)) > eps || M < 4
    error('qam_mod:InvalidOrder', 'M must be a square power of four (4, 16, 64, ...).');
end

bits = bits(:);
if any(bits ~= 0 & bits ~= 1) || mod(numel(bits), k) ~= 0
    error('qam_mod:InvalidBits', 'Input must contain a whole number of binary symbols.');
end
B = reshape(bits, k, []).';
weights = 2.^(m-1:-1:0).';
grayI = B(:,1:m) * weights;
grayQ = B(:,m+1:end) * weights;
indexI = gray_to_binary(grayI);
indexQ = gray_to_binary(grayQ);
levelsI = 2*indexI - (L-1);
levelsQ = 2*indexQ - (L-1);
normalizer = sqrt((2/3) * (M-1));
symbols = (levelsI + 1i*levelsQ) / normalizer;
symbols = symbols(:);
end

function binaryValue = gray_to_binary(grayValue)
value = uint32(grayValue);
shifted = bitshift(value, -1);
while any(shifted(:) ~= 0)
    value = bitxor(value, shifted);
    shifted = bitshift(shifted, -1);
end
binaryValue = double(value);
end
