function [state_dot] = min_dae(state_start,alpha)
    
    state_dot = [0; 0];
        
    B = [ 1 1;
        0 1];
    
    x = state_start(1:2);
    %z = state_start(3:4);
    
    
    A = [ x(1)^2 x(2);
        0 x(2)^2];
    
    A = A + alpha*eye(2);
    epsilon  =  0;
    
    z = inv(A)* x/10;
    state_dot(1:2) =-B * x  - z;  


end
