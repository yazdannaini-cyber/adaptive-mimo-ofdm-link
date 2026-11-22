function decoded = packet_crc_decode_v15_2(receivedEncodedBits,packet)
%PACKET_CRC_DECODE_V15_2 Demap a burst codeword, FEC-decode, and verify CRC.
rx = double(receivedEncodedBits(:));
if ~isstruct(packet) || ~all(isfield(packet,{'CodingMode','PayloadBits','PayloadLength','EncodedBits'}))
    error('packet_crc_decode_v15_2:Packet','Invalid packet metadata.');
end
if numel(rx) ~= numel(packet.EncodedBits) || any(~isfinite(rx)) || ...
        any(rx ~= 0 & rx ~= 1)
    error('packet_crc_decode_v15_2:ReceivedBits','Received codeword length or bit values are invalid.');
end
switch lower(char(packet.CodingMode))
    case 'uncoded-crc'
        protectedHat = rx;
    case 'conv-crc-1/2'
        trellisHat = viterbi_decode_rate12(rx);
        if numel(trellisHat) < 18
            error('packet_crc_decode_v15_2:DecodedLength','Decoded trellis sequence is too short.');
        end
        protectedHat = trellisHat(1:end-2); % strip terminated zero tail
    otherwise
        error('packet_crc_decode_v15_2:CodingMode','Unsupported coding mode in packet metadata.');
end
if numel(protectedHat) ~= packet.PayloadLength+16
    error('packet_crc_decode_v15_2:ProtectedLength','Decoded protected payload has an unexpected length.');
end
payloadHat = protectedHat(1:packet.PayloadLength);
crcValid = crc16_ccitt_verify(protectedHat);
decoded.PayloadBits = payloadHat;
decoded.PayloadErrors = nnz(payloadHat ~= packet.PayloadBits(:));
decoded.PayloadBER = decoded.PayloadErrors/packet.PayloadLength;
decoded.CRCValid = logical(crcValid);
decoded.PacketError = ~decoded.CRCValid;
decoded.UndetectedPayloadErrors = decoded.CRCValid && decoded.PayloadErrors > 0;
decoded.CodedBitErrors = nnz(rx ~= packet.EncodedBits(:));
decoded.CodedBER = decoded.CodedBitErrors/numel(rx);
decoded.ReceivedEncodedBits = numel(rx);
decoded.PayloadLength = packet.PayloadLength;
end
