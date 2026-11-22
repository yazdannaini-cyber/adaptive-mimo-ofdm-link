function [selectedIndex, targetFeasible] = select_qam_by_predicted_ber(predictedBER, targetBER)
%SELECT_QAM_BY_PREDICTED_BER Choose the highest-rate eligible QAM mode.
% predictedBER must be ordered from the lowest to highest modulation order.
% When no mode meets targetBER, index 1 is selected as a minimum-service
% fallback and targetFeasible is false. The fallback is NOT target-compliant.

if isempty(predictedBER) || any(~isfinite(predictedBER(:))) || ...
        any(predictedBER(:) < 0) || any(predictedBER(:) > 0.5) || ...
        ~isscalar(targetBER) || ~isfinite(targetBER) || targetBER <= 0 || targetBER >= 0.5
    error('select_qam_by_predicted_ber:InvalidInput', ...
        'Predicted BER values and target BER must lie in [0,0.5].');
end
eligible = find(predictedBER(:).' <= targetBER);
if isempty(eligible)
    selectedIndex = 1;
    targetFeasible = false;
else
    selectedIndex = eligible(end);
    targetFeasible = true;
end
end
