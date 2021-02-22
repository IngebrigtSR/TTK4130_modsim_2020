%Task 2a) and 2b)
lambda = -2;
t_step = 0.2;
tf = 2;
T = 0:t_step:tf;

f_2b = @(t,x) lambda*x;
J = @(x) -2;

x_0 = 1;

x = ImplicitEulerTemplate(f_2b, J, T, x_0);

%Real Solution
Sol = x_0*exp(lambda*T);

figure(10)
plot(T, Sol, '-', T, x, '-');
legend('Actual Solution', 'Implicit Euler');
title('fig 4: Implicit Euler compared to Actual Solution');