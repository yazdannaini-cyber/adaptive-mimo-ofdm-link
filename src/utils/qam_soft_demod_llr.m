function llr = qam_soft_demod_llr(symbols,M,noiseVariance)
%QAM_SOFT_DEMOD_LLR Max-log bit LLRs for Gray-coded square QAM.
% Positive LLR favors bit 0; negative LLR favors bit 1.
% noiseVariance is E{|n|^2} at the gain-corrected equalizer output and may be
% scalar or one value per input symbol. Output ordering matches qam_demod.
validateattributes(M,{'numeric'},{'scalar','integer','positive'});
k=log2(M); L=sqrt(M); m=log2(L);
if abs(k-round(k))>eps || abs(m-round(m))>eps || M<4
    error('qam_soft_demod_llr:InvalidOrder','M must be a square power of four.');
end
z=symbols(:);
if isempty(z) || any(~isfinite(real(z))) || any(~isfinite(imag(z)))
    error('qam_soft_demod_llr:InvalidSymbols','symbols must be nonempty and finite.');
end
if ~isnumeric(noiseVariance) || ~isreal(noiseVariance) || isempty(noiseVariance) || ...
        any(~isfinite(noiseVariance(:))) || any(noiseVariance(:)<=0) || ...
        ~(isscalar(noiseVariance) || numel(noiseVariance)==numel(z))
    error('qam_soft_demod_llr:InvalidNoise', ...
        'noiseVariance must be positive finite scalar or match the symbol count.');
end
if isscalar(noiseVariance), nv=repmat(double(noiseVariance),numel(z),1);
else, nv=double(noiseVariance(:)); end

% Build the exact bit labels and constellation in qam_mod's ordering.
labels=zeros(M,k);
for q=0:M-1
    for b=1:k
        labels(q+1,b)=double(bitget(uint32(q),k-b+1));
    end
end
constellation=qam_mod(reshape(labels.',[],1),M);
d2=abs(z-constellation.').^2;
llrMatrix=zeros(numel(z),k);
for b=1:k
    d0=min(d2(:,labels(:,b)==0),[],2);
    d1=min(d2(:,labels(:,b)==1),[],2);
    llrMatrix(:,b)=(d1-d0)./max(nv,1e-12);
end
llrMatrix=max(min(llrMatrix,60),-60);
llr=reshape(llrMatrix.',[],1);
end
