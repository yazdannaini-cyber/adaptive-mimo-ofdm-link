function decodedBits = viterbi_decode_rate12_soft(receivedLLR)
%VITERBI_DECODE_RATE12_SOFT Soft-decision Viterbi for terminated K=3 (7,5).
% LLR convention: positive values favor coded bit 0, negative favor bit 1.
% Branch costs are negative-log-likelihood max-log metrics.
r=double(receivedLLR(:));
if isempty(r) || mod(numel(r),2)~=0 || any(~isfinite(r)) || ~isreal(r)
    error('viterbi_decode_rate12_soft:InvalidLLR', ...
        'receivedLLR must be a finite real even-length vector.');
end
nSteps=numel(r)/2;
pathMetric=inf(4,1); pathMetric(1)=0;
predecessor=zeros(nSteps,4,'uint8');
survivorInput=zeros(nSteps,4,'uint8');
for t=1:nSteps
    nextMetric=inf(4,1);
    pair=r(2*t-1:2*t);
    cost0=max(-pair,0); % penalty if hypothesized bit is 0
    cost1=max(pair,0);  % penalty if hypothesized bit is 1
    for state=0:3
        if ~isfinite(pathMetric(state+1)), continue; end
        prev1=bitget(state,2); prev2=bitget(state,1);
        for inputBit=0:1
            reg=4*inputBit+2*prev1+prev2;
            expected=[mod(bitget(reg,3)+bitget(reg,2)+bitget(reg,1),2); ...
                      mod(bitget(reg,3)+bitget(reg,1),2)];
            branch=sum((expected==0).*cost0+(expected==1).*cost1);
            nextState=2*inputBit+prev1;
            metric=pathMetric(state+1)+branch;
            if metric<nextMetric(nextState+1)
                nextMetric(nextState+1)=metric;
                predecessor(t,nextState+1)=uint8(state);
                survivorInput(t,nextState+1)=uint8(inputBit);
            end
        end
    end
    pathMetric=nextMetric;
end
if ~isfinite(pathMetric(1))
    error('viterbi_decode_rate12_soft:NoPath','No terminated trellis path reached state zero.');
end
decodedBits=zeros(nSteps,1); state=0;
for t=nSteps:-1:1
    decodedBits(t)=double(survivorInput(t,state+1));
    state=double(predecessor(t,state+1));
end
end
