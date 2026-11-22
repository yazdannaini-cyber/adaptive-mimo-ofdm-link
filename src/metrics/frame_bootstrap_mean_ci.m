function ci = frame_bootstrap_mean_ci(frameValues, numBootstrap, confidenceLevel, seed)
%FRAME_BOOTSTRAP_MEAN_CI Percentile bootstrap CI for a frame-level mean.
% Each frame-level value is resampled as a whole. No Statistics Toolbox is
% required. The caller should provide one independent value per frame.

frameValues = double(frameValues(:));
if isempty(frameValues) || ~isreal(frameValues) || any(~isfinite(frameValues)) || ...
        ~isscalar(numBootstrap) || ~isfinite(numBootstrap) || ...
        numBootstrap < 200 || numBootstrap ~= floor(numBootstrap) || ...
        ~isscalar(confidenceLevel) || ~isreal(confidenceLevel) || ~isfinite(confidenceLevel) || ...
        confidenceLevel <= 0 || confidenceLevel >= 1 || ...
        ~isscalar(seed) || ~isreal(seed) || ~isfinite(seed) || seed < 0 || seed ~= floor(seed)
    error('frame_bootstrap_mean_ci:InvalidInput', ...
        'Values must be finite and bootstrap settings must be valid.');
end

nFrames = numel(frameValues);
oldRng = rng;
restoreRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(seed,'twister');
indices = randi(nFrames,nFrames,numBootstrap);
bootstrapMeans = sort(mean(frameValues(indices),1));
alpha = 1-confidenceLevel;
ci = [local_percentile(bootstrapMeans,alpha/2), ...
      local_percentile(bootstrapMeans,1-alpha/2)];
end

function value = local_percentile(sortedValues,p)
n = numel(sortedValues);
if n == 1
    value = sortedValues(1);
    return;
end
position = 1+(n-1)*p;
lo = floor(position);
hi = ceil(position);
if lo == hi
    value = sortedValues(lo);
else
    fraction = position-lo;
    value = sortedValues(lo)*(1-fraction)+sortedValues(hi)*fraction;
end
end
