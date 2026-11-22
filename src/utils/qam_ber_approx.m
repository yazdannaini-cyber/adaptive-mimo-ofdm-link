function ber = qam_ber_approx(sinrLinear, M)
%QAM_BER_APPROX Approximate Gray-coded square-QAM BER in AWGN.
% sinrLinear is post-detection symbol SINR (Es/N0). The approximation is
%   Pb ~= (4/log2(M))*(1-1/sqrt(M))*Q(sqrt(3*gamma/(M-1))).
% It is used as a link-adaptation predictor, not as a replacement for the
% Monte-Carlo BER measured from decoded bits.

if ~isscalar(M) || isnan(M) || isinf(M) || M < 4 || abs(log2(M)-round(log2(M))) > eps || ...
        abs(log2(sqrt(M))-round(log2(sqrt(M)))) > eps
    error('qam_ber_approx:InvalidModulation','M must be a square power-of-four QAM order.');
end
if any(isnan(sinrLinear(:))) || any(isinf(sinrLinear(:))) || any(sinrLinear(:) < 0)
    error('qam_ber_approx:InvalidSINR','SINR values must be finite and nonnegative.');
end
k = log2(M);
argument = sqrt((3/(M-1))*sinrLinear);
qValue = 0.5*erfc(argument/sqrt(2));
ber = (4/k)*(1-1/sqrt(M))*qValue;
ber = min(max(ber,0),0.5);
end
