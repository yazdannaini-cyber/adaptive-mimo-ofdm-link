function choice = select_mimo_link_mode(predictedBER, modOrders, targetBER, spectralEfficiency)
%SELECT_MIMO_LINK_MODE Jointly select spatial mode and QAM order.
% Rows of predictedBER/spectralEfficiency are:
%   1: 2x2 spatial multiplexing (two streams)
%   2: 2x2 Alamouti STBC (rate-one code)
% Columns correspond to modOrders in ascending order.
% Among candidates predicted to meet targetBER, choose the highest spectral
% efficiency; ties are resolved by lower predicted BER. If none is feasible,
% choose the candidate with minimum predicted BER and flag targetFeasible=false.

modOrders = double(modOrders(:).');
if isempty(modOrders) || numel(unique(modOrders)) ~= numel(modOrders) || ...
        any(diff(modOrders) <= 0) || any(modOrders < 4) || ...
        any(abs(log2(modOrders)-round(log2(modOrders))) > eps) || ...
        any(abs(log2(sqrt(modOrders))-round(log2(sqrt(modOrders)))) > eps)
    error('select_mimo_link_mode:InvalidModulations', ...
        'modOrders must be strictly increasing square power-of-four QAM orders.');
end
if ~isnumeric(predictedBER) || ~isequal(size(predictedBER),[2 numel(modOrders)]) || ...
        any(~isfinite(predictedBER(:))) || any(predictedBER(:) < 0) || any(predictedBER(:) > 0.5)
    error('select_mimo_link_mode:InvalidPredictedBER', ...
        'predictedBER must be a finite 2-by-numel(modOrders) matrix in [0,0.5].');
end
if ~isnumeric(spectralEfficiency) || ~isequal(size(spectralEfficiency),size(predictedBER)) || ...
        any(~isfinite(spectralEfficiency(:))) || any(spectralEfficiency(:) <= 0)
    error('select_mimo_link_mode:InvalidEfficiency', ...
        'spectralEfficiency must be a positive finite matrix matching predictedBER.');
end
if ~isscalar(targetBER) || ~isfinite(targetBER) || targetBER <= 0 || targetBER >= 0.5
    error('select_mimo_link_mode:InvalidTarget','targetBER must lie strictly between 0 and 0.5.');
end

feasible = predictedBER <= targetBER;
if any(feasible(:))
    candidateEfficiency = spectralEfficiency;
    candidateEfficiency(~feasible) = -Inf;
    bestEfficiency = max(candidateEfficiency(:));
    candidates = find(candidateEfficiency >= bestEfficiency-1e-12);
    [~,bestLocal] = min(predictedBER(candidates));
    linearIndex = candidates(bestLocal);
    targetFeasible = true;
    selectionReason = 'highest-spectral-efficiency-among-predicted-feasible-candidates';
else
    bestBER = min(predictedBER(:));
    candidates = find(predictedBER <= bestBER+1e-12);
    [~,bestLocal] = max(spectralEfficiency(candidates));
    linearIndex = candidates(bestLocal);
    targetFeasible = false;
    selectionReason = 'minimum-predicted-BER-fallback-target-not-met';
end
[modeIndex,modulationIndex] = ind2sub(size(predictedBER),linearIndex);
modeNames = {'spatial-multiplexing','alamouti-stbc'};
choice.ModeIndex = modeIndex;
choice.Mode = modeNames{modeIndex};
choice.ModulationIndex = modulationIndex;
choice.ModulationOrder = modOrders(modulationIndex);
choice.PredictedBER = predictedBER(modeIndex,modulationIndex);
choice.SpectralEfficiency_bpsHz = spectralEfficiency(modeIndex,modulationIndex);
choice.TargetFeasible = targetFeasible;
choice.SelectionReason = selectionReason;
end
