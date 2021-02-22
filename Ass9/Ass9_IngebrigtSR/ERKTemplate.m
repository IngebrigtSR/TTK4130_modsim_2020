function x = ERKTemplate(ButcherArray, f, T, x0)
    % Returns the iterations of an ERK method
    % ButcherArray: Struct with the ERK's Butcher array
    % f: Function handle
    %    Vector field of ODE, i.e., x_dot = f(t,x)
    % T: Vector of time points, 1 x Nt
    % x0: Initial state, Nx x 1
    % x: ERK iterations, Nx x Nt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define variables
    % Allocate space for iterations (x) and k1,k2,...,kNstage
    % It is recommended to allocate a matrix K for all kj, i.e.
    % K = [k1 k2 ... kNstage]
    
    A = ButcherArray.A;    
    c = ButcherArray.c(:);
    b = ButcherArray.b(:);
    
    Nstage = size(A,1);
    Nt = length(T);
    Nx = length(x0);
    
    K = zeros(Nx, Nstage);
    x = zeros(Nx, Nt);
    xt = x0;
    
    dT = diff(T);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x(:,1) = x0;
    % Loop over time points
    for nt=2:Nt
        
        t = T(nt-1);
        dt = dT(nt-1);
        K(:,1) = f(t,xt);
        for nstage=2:Nstage
            a = A(nstage, 1:nstage-1)';
            K(:,nstage) = f(t + dt*c(nstage), xt + dt*(K(:,1:nstage-1)*a));
        end
        xt = xt + dt*(K*b);
        x(:,nt) = xt;
    end
        
        %Mine (the following code) was weird, so used the one from the
        %solution for Ass7 instead
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Update variables
%         x_k = x(:,nt-1);
%         K(:,1) = f(T(nt), x_k);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Loop that calculates k1,k2,...,kNstage
%         for nstage=2:Nstage
%             ksum = 0;
%             for i=1:nstage-1
%                 ksum = ksum + A(nstage,i)*K(:,i);
%             end
%             K(:,nstage) = f(T(nt), x_k+dT*ksum);
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Calculate and save next iteration value x_t
%         xsum = 0;
%         for m=1:Nstage
%            xsum = xsum + b(m)*K(:,m);
%         end
%         x(:,nt) = x_k + dT*xsum;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end