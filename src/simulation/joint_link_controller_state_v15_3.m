function state = joint_link_controller_state_v15_3()
%JOINT_LINK_CONTROLLER_STATE_V15_3 Create empty per-action CRC feedback state.
state.Successes=zeros(2,2,2);
state.Failures=zeros(2,2,2);
state.TotalObservations=zeros(2,2,2);
state.TotalBursts=0;
end
