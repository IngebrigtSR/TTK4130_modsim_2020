clear all;
close all;
clc;

x_start = [1; 1];
y_start = [1; 1];
state_start = [x_start; y_start];
time_finish = 10;
alpha = 0;
epsilon  =  0.000001;

[tsim,xsim]  = ode15s(@(t,state)min_ode(state,alpha, epsilon),[0,time_finish], state_start);

[tsim_dae,xsim_dae]  = ode15s(@(t,state_dae)min_dae(state_dae,alpha),[0,time_finish], x_start);


plot(tsim_dae, xsim_dae(:,1), 'bx')   
hold on
plot(tsim, xsim(:,1))
legend('DAE', 'ODE');

