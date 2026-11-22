function decodedBits = viterbi_decode_rate12(receivedBits)
%VITERBI_DECODE_RATE12 Hard-decision Viterbi decoder for terminated K=3 (7,5) code.
% The terminal state is constrained to zero, matching conv_encode_rate12.
r = double(receivedBits(:));
if isempty(r) || mod(numel(r),2) ~= 0 || any(~isfinite(r)) || ...
        any(r ~= 0 & r ~= 1)
    error('viterbi_decode_rate12:InvalidBits','receivedBits must be an even-length binary vector.');
end
nSteps = numel(r)/2;
pathMetric = inf(4,1); pathMetric(1) = 0;
predecessor = zeros(nSteps,4,'uint8');
survivorInput = zeros(nSteps,4,'uint8');
for t = 1:nSteps
    nextMetric = inf(4,1);
    pair = r(2*t-1:2*t);
    for state = 0:3
        if ~isfinite(pathMetric(state+1)), continue; end
        prev1 = bitget(state,2); prev2 = bitget(state,1);
        for inputBit = 0:1
            reg = 4*inputBit + 2*prev1 + prev2;
            expected = [mod(bitget(reg,3)+bitget(reg,2)+bitget(reg,1),2); ...
                        mod(bitget(reg,3)+bitget(reg,1),2)];
            nextState = 2*inputBit + prev1;
            metric = pathMetric(state+1) + nnz(expected ~= pair);
            if metric < nextMetric(nextState+1)
                nextMetric(nextState+1) = metric;
                predecessor(t,nextState+1) = uint8(state);
                survivorInput(t,nextState+1) = uint8(inputBit);
            end
        end
    end
    pathMetric = nextMetric;
end
if ~isfinite(pathMetric(1))
    error('viterbi_decode_rate12:NoPath','No terminated trellis path reached state zero.');
end
decodedBits = zeros(nSteps,1);
state = 0;
for t = nSteps:-1:1
    decodedBits(t) = double(survivorInput(t,state+1));
    state = double(predecessor(t,state+1));
end
end
