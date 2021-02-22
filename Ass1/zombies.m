%% Initiate all the fucking variables
V.a = 1.4e-6;
V.b = 3.1e-8;
V.b_d = 5.6e-16;
V.d = 2.8e-8;
V.i = 2.6e-6;
V.n = 1.4e-6;
V.r = 2.8e-7;
V.q_i = 2.7e-6;
V.q_z = 2.7e-6;
V.d_q = 2.8e-5;

H0 = (V.b-V.d)/V.b_d;
I0 = 0;
Z0 = 0;
D0 = 0;
Q0 = 0;

tspan = linspace(0,100*24*60*60, 10000);
opts = odeset('RelTol', 0.01, 'AbsTol', 0.01, 'InitialStep', 1);

%% without Q
[t,y] = ode15s(@(t,y) noquarfunc(t,y,V),tspan,[H0 I0 Z0 D0], opts);
figure
plot(t/60/60/24,y);
title("Uten karantene");
xlabel('Dager');
legend('Healthy', 'Infected', 'Zombie', 'Dead');

%% with Q
[t,y] = ode15s(@(t,y) quarfunc(t,y,V), tspan,[H0 I0 Z0 D0 Q0], opts);
figure
plot(t/60/60/24,y);
title("Med karantene");
xlabel('Dager');
legend('Healthy', 'Infected', 'Zombie', 'Dead', 'Quarantined')