function metrics = bit_error_metrics(bitErrors, totalBits, blockErrors, totalBlocks)
%BIT_ERROR_METRICS BER, BLER, and approximate 95% Wilson confidence interval.
if totalBits <= 0 || totalBlocks <= 0
    error('bit_error_metrics:InvalidCounts', 'Totals must be positive.');
end
p = bitErrors / totalBits;
z = 1.96;
denominator = 1 + z^2/totalBits;
center = (p + z^2/(2*totalBits)) / denominator;
halfWidth = z * sqrt(p*(1-p)/totalBits + z^2/(4*totalBits^2)) / denominator;
metrics.BER = p;
metrics.BER_CI95_Lower = max(0, center-halfWidth);
metrics.BER_CI95_Upper = min(1, center+halfWidth);
metrics.BLER = blockErrors / totalBlocks;
metrics.BitErrors = bitErrors;
metrics.TotalBits = totalBits;
metrics.BlockErrors = blockErrors;
metrics.TotalBlocks = totalBlocks;
end
