function [ state_dot ] = SatelliteDynamics( t, x, parameters )

    % Code your equations here...
    
    
    
   G = 6.67*10^(-11);
    m_T = 5.97*10^(24);
    r = 6356*10^3 + 36*10^6;
    
%   x = [pos;
%          R;
%          vel;
%          ang_vel];

    pos = x(1:3);
    
    R = reshape(x(4:12),3,3);
    
    vel = x(13:15);
    
    ang_vel = x(16:18);
    
    
%   parameters = [M;
%                 m_T;]

    omega_skew =  [0 -ang_vel(3) ang_vel(2) ; ang_vel(3) 0 -ang_vel(1) ; -ang_vel(2) ang_vel(1) 0 ]; 
    r_dot = vel;
    R_dot = R * omega_skew;
    acc = -(G*m_T*r)*pos(1:3)/(r^2);
    ang_acc = [0; 0; 0];%- (parameters(1) * eye(3))' * omega_skew * parameters(1)* eye(3) * ang_vel;
    
%      state_dot(1:3) = state(13:15);
%     R_dot = reshape(state(4:12),3,3)*omega_ab_skew;
%     state_dot(4:12) = reshape(R_dot, 9, 1);
%     state_dot(13:15) = -(G*m_T*r)*state(1:3)/(r^2);
%     state_dot(16:18) = [0; 0; 0];
%     state_dot = state_dot';   
    % The code must return in the order you selected, e.g.:
    state_dot = [ r_dot;
                    reshape(R_dot,9,1);
                    acc;
                    ang_acc];
    
    %    state_dot =  [velocity;
    %                  orientation_dot;
    %                  acceleration (ac);
    %                  angular acceleration (omega dot)];

end
