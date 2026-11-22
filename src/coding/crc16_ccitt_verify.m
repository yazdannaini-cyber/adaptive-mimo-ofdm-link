function isValid = crc16_ccitt_verify(codeword)
%CRC16_CCITT_VERIFY Verify a CRC-16/CCITT-FALSE appended binary codeword.
codeword = double(codeword(:));
if numel(codeword) < 17 || any(~isfinite(codeword)) || ...
        any(codeword ~= 0 & codeword ~= 1)
    error('crc16_ccitt_verify:InvalidBits','codeword must contain payload and 16 binary CRC bits.');
end
expected = crc16_ccitt_append(codeword(1:end-16));
isValid = isequal(expected,codeword);
end
