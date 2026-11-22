function decoded = packet_soft_decode_v15_2(receivedLLR,packet)
%PACKET_SOFT_DECODE_V15_2 Soft demap, optional soft Viterbi, then CRC check.
llr=double(receivedLLR(:));
if ~isstruct(packet) || ~all(isfield(packet,{'CodingMode','PayloadBits','PayloadLength','EncodedBits'}))
    error('packet_soft_decode_v15_2:Packet','Invalid packet metadata.');
end
if numel(llr)~=numel(packet.EncodedBits) || ~isreal(llr) || any(~isfinite(llr))
    error('packet_soft_decode_v15_2:LLR','LLR length must match encoded packet and be finite real values.');
end
hardBits=double(llr<0);
switch lower(char(packet.CodingMode))
    case 'uncoded-crc'
        protectedHat=hardBits;
    case 'conv-crc-1/2'
        trellisHat=viterbi_decode_rate12_soft(llr);
        if numel(trellisHat)<18
            error('packet_soft_decode_v15_2:DecodedLength','Decoded trellis sequence is too short.');
        end
        protectedHat=trellisHat(1:end-2);
    otherwise
        error('packet_soft_decode_v15_2:CodingMode','Unsupported packet coding mode.');
end
if numel(protectedHat)~=packet.PayloadLength+16
    error('packet_soft_decode_v15_2:ProtectedLength','Decoded protected payload has an unexpected length.');
end
payloadHat=protectedHat(1:packet.PayloadLength);
crcValid=crc16_ccitt_verify(protectedHat);
decoded.PayloadBits=payloadHat;
decoded.PayloadErrors=nnz(payloadHat~=packet.PayloadBits(:));
decoded.PayloadBER=decoded.PayloadErrors/packet.PayloadLength;
decoded.CRCValid=logical(crcValid);
decoded.PacketError=~decoded.CRCValid;
decoded.UndetectedPayloadErrors=decoded.CRCValid && decoded.PayloadErrors>0;
decoded.CodedBitErrors=nnz(hardBits~=packet.EncodedBits(:));
decoded.CodedBER=decoded.CodedBitErrors/numel(hardBits);
decoded.ReceivedEncodedBits=numel(hardBits);
decoded.PayloadLength=packet.PayloadLength;
end
