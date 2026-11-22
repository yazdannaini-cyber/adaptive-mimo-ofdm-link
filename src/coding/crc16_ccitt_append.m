function codeword = crc16_ccitt_append(payloadBits)
%CRC16_CCITT_APPEND Append CRC-16/CCITT-FALSE (poly 0x1021, init 0xFFFF).
% Input and output bits are column vectors, MSB-first within each byte/CRC word.
payloadBits = double(payloadBits(:));
if isempty(payloadBits) || any(~isfinite(payloadBits)) || ...
        any(payloadBits ~= 0 & payloadBits ~= 1)
    error('crc16_ccitt_append:InvalidBits','payloadBits must be a nonempty binary vector.');
end
crc = uint16(hex2dec('FFFF'));
poly = uint16(hex2dec('1021'));
for k = 1:numel(payloadBits)
    feedback = xor(bitget(crc,16) ~= 0, payloadBits(k) ~= 0);
    crc = bitshift(crc,1);
    if feedback
        crc = bitxor(crc,poly);
    end
end
crcBits = zeros(16,1);
for k = 1:16
    crcBits(k) = double(bitget(crc,17-k));
end
codeword = [payloadBits;crcBits];
end
