%IRKTemplate(ButcherArray, f, dfdx, T, x0)
clc;
close all;
clear all;
A_IRK = [(1/4)  (1/4-(sqrt(3)/6)); 
     (1/4 + (sqrt(3)/6)) (1/4)];
c_IRK = [(1/2 - (sqrt(3)/6)); 
     (1/2 + (sqrt(3)/6))];
b_IRK = [(1/2) (1/2)];
ButcherExample_IRK = struct('A',A_IRK,'b',b_IRK,'c',c_IRK);

lambda = -2;
t_final = 2;
t_step = 0.4;
T_linspace = linspace(0,t_final,100);
x0 = 1;

f = @(t,x) lambda*x;
J = @(t,x) lambda;
S_anal = exp(lambda*T_linspace);

T = 0:t_step:t_final;

S_IRK = IRKTemplate(ButcherExample_IRK, f, J, T, x0);

% plot(T, S_IRK)
% hold on

A_ERK = [0 0 0 0;
        1/2 0 0 0;
        0 1/2 0 0;
        0 0 1 0];
c_ERK = [0;
         1/2;
         1/2;
         1];
b_ERK = [1/6;
         1/3;
         1/3;
         1/6];

ButcherExample_ERK = struct('A',A_ERK,'b',b_ERK,'c',c_ERK);

S_ERK = ERKTemplate(ButcherExample_ERK, f, T, x0);

% plot(T, S_ERK, T_linspace, S_anal, '--');
% legend('IRK', 'ERK', 'Actual/Analytic');
% title('Figure 3')

%%Task 2)

x_0 = 2;
xdot_0 = 0;
x_init = [x_0;
          xdot_0];

dt_k = 0.01;
tf = 10;
T2 = 0:dt_k:tf;

x_d = 1.32;
kappa = 2.4;
g = 9.81;
m = 200;

f2 = @(t,x) [x(2);
            (-g*(1-(x_d/x(1))^kappa))];
            
J2 = @(t,x) [0 1;
            (-(g*kappa*x_d*(x_d/x(1))^(kappa-1))/x(1)^2) 0];
        
A_EE = 0;
c_EE = 0;
b_EE = 1;
ButcherExample_EE = struct('A',A_EE,'b',b_EE,'c',c_EE);
            
A_G2 = 1/2;
c_G2 = 1/2;
b_G2 = 1;

ButcherExample_G2 = struct('A',A_G2,'b',b_G2,'c',c_G2);

S2_IRK = IRKTemplate(ButcherExample_IRK, f2, J2, T2, x_init)';
S2_EE = ERKTemplate(ButcherExample_EE, f2, T2, x_init)';
S2_G2 = IRKTemplate(ButcherExample_G2, f2, J2, T2, x_init)';

% subplot(211)
% plot(T2, S2_IRK(:,1))
% hold on
% plot(T2, S2_EE(:,1));
% hold on
% plot(T2, S2_G2(:,1));
% hold on
% legend('IRK', 'EE', 'G2');
% title('Figure 4');
% 
% subplot(212)
% plot(T2, S2_IRK(:,2))
% hold on
% plot(T2, S2_EE(:,2));
% hold on
% plot(T2, S2_G2(:,2));
% hold on
% legend('IRK', 'EE', 'G2');
% title('Figure 5');

%%Task 3)
p  = sym('p', [3,1]);
v  = sym('v', [3,1]);
  
t = sym('t');
z = sym('z');

dp = sym('dp', [3,1]);
dv = sym('dv', [3,1]);

X  = [p; 
      v];
dX = [dp; 
      dv];

X_0 = [1 0 0 0 0 1]'; %#2 should be 0
Z_0 = 1;
L = 1;
M = 10;

%Cq = 1/2 * (X(1:3)' * X(1:3) - L^2); for 3b)

f = [X(4:6) - dX(1:3);
    -M*g*[0 0 1]' - z*X(1:3) - M*dX(4:6);
    X(1:3)' * dX(4:6) + X(4:6)'*X(4:6)];
    
    %X(1:3)' * dX(4:6) + X(4:6)'*X(4:6)];
    



j_dX = jacobian(f, dX);
j_X = jacobian(f, X);
j_Z = jacobian(f, z);


J_dX = matlabFunction(j_dX, 'Vars', {[dX],[X],z,t});
J_X = matlabFunction(j_X, 'Vars', {[dX],[X],z,t});
J_Z = matlabFunction(j_Z, 'Vars', {[dX],[X],z,t});
f = matlabFunction(f, 'Vars', {[dX],[X],z,t});


dT3 = 0.01;
T3 = [0:dT3:30];



[X, dX, z] = RKDAE(ButcherExample_IRK, f, J_dX, J_X, J_Z, T3, X_0, Z_0);


% j_Cq = jacobian(Cq,X(1:3,:)');
% Cq = matlabFunction(Cq, 'Vars', {[dX],[X],z,t});
% J_Cq = matlabFunction(j_Cq, 'Vars', {[dX],[X],z,t});
% Sol_Cq = IRKTemplate(ButcherExample_IRK, Cq, j_Cq, T3, X_0);

% plot(T3,z);

% subplot(211)
% plot(T3, X(1:3,:)');
% legend('X','Y','Z');
% title('Figure 6: 3D-position of the pendulum')
% 
% subplot(212)
% plot(T3, X(4:6,:)');
% legend('dX','dY','dZ');
% title('Figure 7: Component speeds of the pendulum')

