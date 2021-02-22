function [ state_dot ] = Kinematics( t, state, parameters )
    % state_dot is time derivative of your state.
    %for 2a)
    %state_dot = inv(M)*parameters 
    
    %for 2b)
    omega_skew = [0 -parameters(3) parameters(2); 
                    parameters(3) 0 -parameters(1); 
                    -parameters(2) parameters(1) 0];
    dRba = reshape(state, [3,3])*omega_skew;
    state_dot = reshape(dRba, 9, 1);
    % Hints:
    % - "parameters" allows you to pass some parameters to the "Kinematic" function.
    % - "state" will contain representations of the solid orientation (SO(3)).
    % - use the "reshape" function to turn a matrix into a vector or vice-versa.

    % Code your equations here...
end
