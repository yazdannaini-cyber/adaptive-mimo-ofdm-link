function [choice,diagnostics] = select_joint_mimo_mcs_fec_v15_3( ...
    prediction,nDataSymbols,nFFT,cpLength,targetPER,state,options)
%SELECT_JOINT_MIMO_MCS_FEC_V15_3 Joint mode/QAM/FEC selector with CRC feedback.
% Eight actions are evaluated: SM/Alamouti x QPSK/16-QAM x uncoded/Conv(7,5).
% A conservative link-model PER prior is updated with causal CRC ACK/NACK
% observations from previously completed bursts. The decision maximizes the
% posterior expected CRC-accepted payload goodput subject to a one-sided
% 95% PER constraint. If no action is feasible, maximum expected goodput is
% selected and the feasibility flag remains false.
%
% The PER prior is an online link-adaptation proxy, not a hard guarantee.
% It is corrected by packet-level feedback during the simulation.

if nargin < 6 || isempty(state)
    state = joint_link_controller_state_v15_3();
end
if nargin < 7 || isempty(options)
    options = struct();
end
if ~isstruct(options)
    error('select_joint_mimo_mcs_fec_v15_3:Options','options must be a struct.');
end
if ~isstruct(prediction) || ~all(isfield(prediction, ...
        {'predictedBER','modulations','postSinrSpatialMultiplexing','postSinrAlamouti'}))
    error('select_joint_mimo_mcs_fec_v15_3:Prediction', ...
        'prediction must contain predictedBER, modulations and post-SINR arrays.');
end
if ~isscalar(nDataSymbols) || ~isfinite(nDataSymbols) || nDataSymbols < 2 || ...
        nDataSymbols ~= floor(nDataSymbols) || mod(nDataSymbols,2) ~= 0 || ...
        ~isscalar(nFFT) || ~isfinite(nFFT) || nFFT < 2 || nFFT ~= floor(nFFT) || ...
        ~isscalar(cpLength) || ~isfinite(cpLength) || cpLength < 0 || ...
        cpLength ~= floor(cpLength) || cpLength >= nFFT
    error('select_joint_mimo_mcs_fec_v15_3:Dimensions', ...
        'Invalid OFDM dimensions or data-symbol count.');
end
if ~isscalar(targetPER) || ~isfinite(targetPER) || targetPER <= 0 || targetPER >= 1
    error('select_joint_mimo_mcs_fec_v15_3:TargetPER', ...
        'targetPER must lie strictly between zero and one.');
end
if ~isstruct(state) || ~all(isfield(state,{'Successes','Failures','TotalObservations'})) || ...
        ~isequal(size(state.Successes),[2 2 2]) || ...
        ~isequal(size(state.Failures),[2 2 2]) || ...
        ~isfield(state,'TotalObservations') || ~isequal(size(state.TotalObservations),[2 2 2]) || ...
        any(~isfinite(state.TotalObservations(:))) || any(state.TotalObservations(:)<0) || ...
        any(state.TotalObservations(:)~=floor(state.TotalObservations(:))) || ...
        any(~isfinite(state.Successes(:))) || any(~isfinite(state.Failures(:))) || ...
        any(state.Successes(:)<0) || any(state.Failures(:)<0) || ...
        any(state.Successes(:)~=floor(state.Successes(:))) || ...
        any(state.Failures(:)~=floor(state.Failures(:))) || ...
        any(state.TotalObservations(:) ~= state.Successes(:)+state.Failures(:))
    error('select_joint_mimo_mcs_fec_v15_3:ControllerState', ...
        'Controller state must contain nonnegative 2-by-2-by-2 success/failure counts.');
end

priorWeight = get_option(options,'priorWeight',4);
fecCodingGainDb = get_option(options,'fecCodingGainDb',2.0);
highBERNoGainThreshold = get_option(options,'highBERNoGainThreshold',0.12);
if ~isscalar(priorWeight)||~isfinite(priorWeight)||priorWeight<=0 || ...
        ~isscalar(fecCodingGainDb)||~isfinite(fecCodingGainDb)||fecCodingGainDb<0 || ...
        ~isscalar(highBERNoGainThreshold)||~isfinite(highBERNoGainThreshold)|| ...
        highBERNoGainThreshold<=0||highBERNoGainThreshold>=0.5
    error('select_joint_mimo_mcs_fec_v15_3:Options', ...
        'Invalid priorWeight, fecCodingGainDb or highBERNoGainThreshold.');
end

mods = double(prediction.modulations(:).');
if ~isnumeric(prediction.postSinrSpatialMultiplexing) || isempty(prediction.postSinrSpatialMultiplexing) || ...
        ~isnumeric(prediction.postSinrAlamouti) || isempty(prediction.postSinrAlamouti) || ...
        any(~isfinite(prediction.postSinrSpatialMultiplexing(:))) || ...
        any(~isfinite(prediction.postSinrAlamouti(:))) || ...
        any(prediction.postSinrSpatialMultiplexing(:)<0) || any(prediction.postSinrAlamouti(:)<0)
    error('select_joint_mimo_mcs_fec_v15_3:PostSINR','Post-equalization SINR arrays must be nonempty, finite and nonnegative.');
end
if size(prediction.postSinrSpatialMultiplexing,1)~=nFFT || ...
        size(prediction.postSinrSpatialMultiplexing,2)~=2 || ...
        size(prediction.postSinrSpatialMultiplexing,3)~=nDataSymbols || ...
        size(prediction.postSinrAlamouti,1)~=nFFT || ...
        size(prediction.postSinrAlamouti,2)~=2 || ...
        size(prediction.postSinrAlamouti,3)~=nDataSymbols/2
    error('select_joint_mimo_mcs_fec_v15_3:PostSINRDimensions', ...
        'Post-SINR arrays must be [nFFT x 2 x nDataSymbols] and [nFFT x 2 x nDataSymbols/2].');
end
if ~isequal(size(prediction.predictedBER),[2 numel(mods)]) || ...
        any(~isfinite(prediction.predictedBER(:))) || ...
        any(prediction.predictedBER(:)<0) || any(prediction.predictedBER(:)>0.5)
    error('select_joint_mimo_mcs_fec_v15_3:PredictedBER', ...
        'predictedBER must be a finite 2-by-numel(modulations) matrix in [0,0.5].');
end
if ~isequal(mods,[4 16])
    error('select_joint_mimo_mcs_fec_v15_3:Modulations', ...
        'V15.3 joint controller currently supports exactly QPSK and 16-QAM.');
end

nActions = 2*numel(mods)*2;
actions = repmat(struct('Mode','','ModulationOrder',NaN,'CodingMode','', ...
    'PredictedPayloadBER',NaN,'ModelPER',NaN,'PosteriorPER',NaN, ...
    'PosteriorPERUpper95',NaN,'PayloadSE_bpsHz',NaN,'ExpectedCRCGoodput_bpsHz',NaN, ...
    'PayloadBits',NaN,'CodeRate',NaN,'CodingGainDb',NaN, ...
    'SuccessCount',0,'FailureCount',0,'TargetFeasible',false, ...
    'ModeIndex',0,'ModulationIndex',0,'CodingIndex',0),nActions,1);

idx = 0;
for modeIndex = 1:2
    if modeIndex==1, modeName='spatial-multiplexing'; else, modeName='alamouti-stbc'; end
    for modIndex = 1:numel(mods)
        M = mods(modIndex); bitsPerSymbol = log2(M);
        if modeIndex==1
            encodedCapacity = nFFT*bitsPerSymbol*2*nDataSymbols;
            sinrValues = prediction.postSinrSpatialMultiplexing(:);
        else
            encodedCapacity = nFFT*bitsPerSymbol*nDataSymbols;
            sinrValues = prediction.postSinrAlamouti(:);
        end
        if mod(encodedCapacity,2)~=0
            error('select_joint_mimo_mcs_fec_v15_3:OddCapacity', ...
                'All candidate encoded capacities must be even.');
        end
        rawBER = prediction.predictedBER(modeIndex,modIndex);
        for codingIndex = 1:2
            idx=idx+1;
            if codingIndex==1
                codingMode='uncoded-crc';
                payloadBits=encodedCapacity-16;
                codeRate=1;
                candidateBER=rawBER;
                codingGain=0;
            else
                codingMode='conv-crc-1/2';
                payloadBits=encodedCapacity/2-18;
                codeRate=0.5;
                if rawBER>=highBERNoGainThreshold
                    codingGain=0;
                elseif rawBER>=0.05
                    codingGain=0.75*fecCodingGainDb;
                else
                    codingGain=fecCodingGainDb;
                end
                if codingGain > 0
                    gainLinear=10^(codingGain/10);
                    candidateBER=mean(qam_ber_approx(sinrValues*gainLinear,M));
                else
                    % Do not claim any FEC benefit outside the supported
                    % coding-gain region: retain the uncoded link BER prior.
                    candidateBER=rawBER;
                end
            end
            if payloadBits<1
                error('select_joint_mimo_mcs_fec_v15_3:PacketCapacity', ...
                    'Candidate payload capacity must be positive.');
            end
            % Packet-error prior assumes independent bit errors only as an
            % initial conservative model. Real CRC outcomes replace it via
            % feedback; it is not reported as a measured PER.
            modelPER = -expm1(payloadBits*log(max(1-candidateBER,realmin)));
            modelPER = min(max(modelPER,0),1);
            successCount=state.Successes(modeIndex,modIndex,codingIndex);
            failureCount=state.Failures(modeIndex,modIndex,codingIndex);
            nObs=successCount+failureCount;
            posteriorPER=(priorWeight*modelPER+failureCount)/(priorWeight+nObs);
            upper95=wilson_upper_fractional(failureCount+priorWeight*modelPER,priorWeight+nObs,1.64485362695147);
            payloadSE=payloadBits/(nFFT*nDataSymbols)*nFFT/(nFFT+cpLength);
            expectedGoodput=payloadSE*(1-posteriorPER);
            actions(idx).Mode=modeName;
            actions(idx).ModulationOrder=M;
            actions(idx).CodingMode=codingMode;
            actions(idx).PredictedPayloadBER=candidateBER;
            actions(idx).ModelPER=modelPER;
            actions(idx).PosteriorPER=posteriorPER;
            actions(idx).PosteriorPERUpper95=upper95;
            actions(idx).PayloadSE_bpsHz=payloadSE;
            actions(idx).ExpectedCRCGoodput_bpsHz=expectedGoodput;
            actions(idx).PayloadBits=payloadBits;
            actions(idx).CodeRate=codeRate;
            actions(idx).CodingGainDb=codingGain;
            actions(idx).SuccessCount=successCount;
            actions(idx).FailureCount=failureCount;
            actions(idx).TargetFeasible=upper95<=targetPER;
            actions(idx).ModeIndex=modeIndex;
            actions(idx).ModulationIndex=modIndex;
            actions(idx).CodingIndex=codingIndex;
        end
    end
end

feasible = find([actions.TargetFeasible]);
if ~isempty(feasible)
    utility = [actions(feasible).ExpectedCRCGoodput_bpsHz];
    bestUtility=max(utility);
    tied=feasible(utility>=bestUtility-1e-12);
    [~,local]=min([actions(tied).PosteriorPERUpper95]);
    selected=tied(local);
    reason='maximum-posterior-CRC-goodput-among-PER-feasible-actions';
else
    utility=[actions.ExpectedCRCGoodput_bpsHz];
    bestUtility=max(utility);
    tied=find(utility>=bestUtility-1e-12);
    [~,local]=min([actions(tied).PosteriorPER]);
    selected=tied(local);
    reason='no-action-meets-PER-target-maximum-expected-goodput-fallback';
end

choice=actions(selected);
choice.ActionIndex=selected;
choice.TargetPER=targetPER;
choice.TargetPERFeasible=choice.TargetFeasible;
choice.SelectionReason=reason;
choice.PredictionSource='SINR-link-abstraction-plus-causal-CRC-feedback';
diagnostics.Actions=actions;
diagnostics.SelectedActionIndex=selected;
diagnostics.FeasibleActionCount=numel(feasible);
diagnostics.TotalFeedbackObservations=sum(state.TotalObservations(:));
diagnostics.Options.PriorWeight=priorWeight;
diagnostics.Options.FECCodingGainDb=fecCodingGainDb;
diagnostics.Options.HighBERNoGainThreshold=highBERNoGainThreshold;
end

function value=get_option(options,name,defaultValue)
if isfield(options,name), value=options.(name); else, value=defaultValue; end
end

function upper=wilson_upper_fractional(errors,total,z)
if total<=0, upper=1; return; end
p=min(max(errors/total,0),1);
den=1+z^2/total;
center=(p+z^2/(2*total))/den;
half=z*sqrt(max(p*(1-p)/total+z^2/(4*total^2),0))/den;
upper=min(1,center+half);
end
