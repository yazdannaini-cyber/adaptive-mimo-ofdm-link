function state = update_joint_link_controller_state_v15_3(state,choice,packetError)
%UPDATE_JOINT_LINK_CONTROLLER_STATE_V15_3 Add one causal packet CRC outcome.
if ~isstruct(state)||~all(isfield(state,{'Successes','Failures','TotalObservations'}))|| ...
        ~isequal(size(state.Successes),[2 2 2])||~isequal(size(state.Failures),[2 2 2])|| ...
        ~isequal(size(state.TotalObservations),[2 2 2])
    error('update_joint_link_controller_state_v15_3:State','Invalid controller state.');
end
if ~isstruct(choice)||~all(isfield(choice,{'Mode','ModulationOrder','CodingMode'}))|| ...
        ~isscalar(packetError)||~(islogical(packetError)||isnumeric(packetError))|| ...
        ~isfinite(packetError)||~ismember(double(packetError),[0 1])
    error('update_joint_link_controller_state_v15_3:Input','Invalid action or packetError flag.');
end
if strcmpi(choice.Mode,'spatial-multiplexing'), modeIndex=1;
elseif strcmpi(choice.Mode,'alamouti-stbc'), modeIndex=2;
else, error('update_joint_link_controller_state_v15_3:Mode','Unknown MIMO mode.'); end
if choice.ModulationOrder==4, modIndex=1;
elseif choice.ModulationOrder==16, modIndex=2;
else, error('update_joint_link_controller_state_v15_3:Modulation','Unsupported modulation.'); end
if strcmpi(choice.CodingMode,'uncoded-crc'), codingIndex=1;
elseif strcmpi(choice.CodingMode,'conv-crc-1/2'), codingIndex=2;
else, error('update_joint_link_controller_state_v15_3:Coding','Unknown coding mode.'); end
if logical(packetError), state.Failures(modeIndex,modIndex,codingIndex)=state.Failures(modeIndex,modIndex,codingIndex)+1;
else, state.Successes(modeIndex,modIndex,codingIndex)=state.Successes(modeIndex,modIndex,codingIndex)+1; end
state.TotalObservations(modeIndex,modIndex,codingIndex)= ...
    state.TotalObservations(modeIndex,modIndex,codingIndex)+1;
if ~isfield(state,'TotalBursts'), state.TotalBursts=0; end
state.TotalBursts=state.TotalBursts+1;
end
