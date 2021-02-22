function [state_dot] = min_ode(state_start, alpha,epsilon)

    state_dot = [0; 0];        
    B = [ 1 1;
        0 1];
    x = state_start(1:2);
    z = state_start(3:4);
    
    A = [ x(1)^2 x(2);
        0 x(2)^2];
    
    A = A + alpha*eye(2);
    

    state_dot(3:4) = (1/10 * x -  A * z) / epsilon;
    
    state_dot(1:2) =-B * x  - z;  
    
    
    

end
