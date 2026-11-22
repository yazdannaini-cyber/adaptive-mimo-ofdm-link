function packet = packet_crc_encode_v15_2(encodedCapacityBits,codingMode,seed)
%PACKET_CRC_ENCODE_V15_2 Build one burst packet with CRC, optionally rate-1/2 FEC.
% codingMode: 'uncoded-crc' or 'conv-crc-1/2'. Output encodedBits always fills
% exactly encodedCapacityBits available QAM bit positions.
if ~isscalar(encodedCapacityBits) || ~isfinite(encodedCapacityBits) || ...
        encodedCapacityBits < 64 || encodedCapacityBits ~= floor(encodedCapacityBits)
    error('packet_crc_encode_v15_2:Capacity','encodedCapacityBits must be an integer >=64.');
end
if ~isscalar(seed) || ~isfinite(seed) || seed < 0 || seed ~= floor(seed)
    error('packet_crc_encode_v15_2:Seed','seed must be a nonnegative integer.');
end
mode = lower(char(codingMode));
switch mode
    case 'uncoded-crc'
        payloadLength = encodedCapacityBits-16;
    case 'conv-crc-1/2'
        if mod(encodedCapacityBits,2) ~= 0
            error('packet_crc_encode_v15_2:OddCapacity','Rate-1/2 code requires an even encoded capacity.');
        end
        payloadLength = encodedCapacityBits/2-16-2; % CRC plus two trellis-tail bits
    otherwise
        error('packet_crc_encode_v15_2:CodingMode','Unsupported coding mode.');
end
if payloadLength < 1
    error('packet_crc_encode_v15_2:InsufficientCapacity','Packet capacity is too small for CRC/FEC overhead.');
end
oldRng = rng; restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(seed,'twister');
payloadBits = randi([0 1],payloadLength,1);
protectedBits = crc16_ccitt_append(payloadBits);
if strcmp(mode,'uncoded-crc')
    encodedBits = protectedBits;
else
    encodedBits = conv_encode_rate12(protectedBits);
end
if numel(encodedBits) ~= encodedCapacityBits
    error('packet_crc_encode_v15_2:InternalLength','Encoded packet does not fill its allocated capacity.');
end
packet.CodingMode = mode;
packet.PayloadBits = payloadBits;
packet.PayloadLength = payloadLength;
packet.ProtectedBits = protectedBits;
packet.EncodedBits = encodedBits;
packet.EncodedCapacityBits = encodedCapacityBits;
packet.CRCBits = 16;
if strcmp(mode,'uncoded-crc'), packet.CodeRate = 1; else, packet.CodeRate = 0.5; end
end
