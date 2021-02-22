clear
clc
% Mass-damper-spring parameters
m = 1;
d = 0.1;
k = 1;
A = [ 0      1;
    -k/m -d/m];
% Mass-damper-spring vector field
fMassDamperSpring = @(t,x) A*x;

% Explicit Euler
A1 = 0;
c1 = 0;
b1 = 1;
RK1 = struct('A',A1,'b',b1,'c',c1);

% RK2
A2 = [0  0; 
     1/2 0];
c2 = [0; 
     1/2];
b2 = [0; 
      1];
RK2 = struct('A',A2,'b',b2,'c',c2);

% RK4
A4 =  [0 0 0 0; 
     1/2 0 0 0; 
     0 1/2 0 0; 
      0 0 1 0];
c4 = [0; 
    1/2; 
    1/2; 
    1];
b4 = [1/6; 
    1/3; 
    1/3; 
    1/6];
RK4 = struct('A',A4,'b',b4,'c',c4);


% Task 1-2 parameters
lambda = -2;
dT = 0.11;
T = 0:dT:25;
x0 = 1;
func = @(t,x) lambda*x;
actual_solution = exp(lambda*T);

% Simulate
X1 = ERKTemplate(RK1,func,T,dT,x0);
X2 = ERKTemplate(RK2,func,T,dT,x0);
X4 = ERKTemplate(RK4,func,T,dT,x0);
X_mass = ERKTemplate(RK2,fMassDamperSpring,T,dT,[1;1]);

% Task 2 plots

%mass damper spring
figure(1)
plot(T,X_mass(1,:),T,X_mass(2,:), '--'); 
legend('Position [m]','Velocity [m/s]');
xlabel('T')
grid on

%RK1,2,4
figure(2)
subplot(3,1,1)
plot(T,X1,T,actual_solution, '--'); 
legend('Explicit Euler','The actual solution');
ylabel('x(t)'); 
xlabel('T'); 
title('Explicit Euler');
grid on

subplot(3,1,2)
plot(T,X2,T,actual_solution, '--'); 
legend('RK2','The actual solution');
ylabel('x(t)'); 
xlabel('T'); 
title('RK2');
grid on

subplot(3,1,3)
plot(T,X4,T,actual_solution, '--'); 
legend('RK4','The actual solution');
ylabel('x(t)'); 
xlabel('T'); 
title('RK4');
grid on

% Task 3 vanderpol
u = 5;
state0 = [2;
          0];
t_final = 25;

[time,statetraj] = ode45(@(t,x)vanderpol(t, x, u),[0 t_final], state0);

vanderpol_func = @(t,x) vanderpol(t, x, u);

x_vdp = ERKTemplate(RK4,vanderpol_func,T,dT,state0);

%Task 3 Plot
figure(3)
subplot(2,1,1)
plot(time,statetraj(:,1),T,x_vdp(1,:), '--');
ylabel('x(t)'); legend('ODE45', 'RK4');
grid on

subplot(2,1,2)
plot(time,statetraj(:,2),T,x_vdp(2,:), '--');
ylabel('y(t)'); xlabel('T'); legend('ODE45', 'RK4');
grid on