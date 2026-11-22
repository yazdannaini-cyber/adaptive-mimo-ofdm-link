function ci = cluster_bootstrap_ber_ci(frameErrors, frameBits, numBootstrap, confidenceLevel, seed)
%CLUSTER_BOOTSTRAP_BER_CI Frame-cluster bootstrap interval for aggregate BER.
% Each frame is resampled as a whole. This preserves the dependence among
% bits that share the same channel realization within a simulated frame.
% The point estimate remains sum(frameErrors)/sum(frameBits); this function
% estimates uncertainty around that ratio using a percentile bootstrap.
%
% Note: if every observed frame has zero errors, the nonparametric bootstrap
% interval is degenerate at zero. Report it alongside the bit-IID Wilson
% interval (which has a nonzero upper bound for zero observed errors).

frameErrors = double(frameErrors(:));
frameBits = double(frameBits(:));
if isempty(frameErrors) || numel(frameErrors) ~= numel(frameBits) || ...
        any(~isfinite(frameErrors)) || any(~isfinite(frameBits)) || ...
        any(frameErrors < 0) || any(frameBits <= 0) || ...
        any(frameErrors ~= floor(frameErrors)) || any(frameBits ~= floor(frameBits)) || ...
        any(frameErrors > frameBits) || ...
        ~isscalar(numBootstrap) || ~isfinite(numBootstrap) || numBootstrap < 200 || ...
        numBootstrap ~= floor(numBootstrap) || ...
        ~isscalar(confidenceLevel) || ~isfinite(confidenceLevel) || ...
        confidenceLevel <= 0 || confidenceLevel >= 1 || ...
        ~isscalar(seed) || ~isfinite(seed) || seed < 0 || seed ~= floor(seed)
    error('cluster_bootstrap_ber_ci:InvalidInput','Invalid frame counts or bootstrap settings.');
end

pointEstimate = sum(frameErrors) / sum(frameBits);
numFrames = numel(frameErrors);
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(seed,'twister');
indices = randi(numFrames,numFrames,numBootstrap);
bootErrors = sum(frameErrors(indices),1);
bootBits = sum(frameBits(indices),1);
bootBER = sort(bootErrors ./ bootBits);
alpha = 1-confidenceLevel;
ci = [local_percentile(bootBER,alpha/2), local_percentile(bootBER,1-alpha/2)];
ci(1) = max(0,min(ci(1),pointEstimate));
ci(2) = min(1,max(ci(2),pointEstimate));
end

function value = local_percentile(sortedValues,p)
% Linear interpolation percentile without Statistics Toolbox dependencies.
n = numel(sortedValues);
if n == 1
    value = sortedValues(1);
    return;
end
position = 1 + (n-1)*p;
lo = floor(position);
hi = ceil(position);
if lo == hi
    value = sortedValues(lo);
else
    fraction = position-lo;
    value = sortedValues(lo)*(1-fraction) + sortedValues(hi)*fraction;
end
end
