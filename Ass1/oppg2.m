
ode1 = @(t,x) x.^2;
t1 = [0 5];
x1_0 = 1;
ode2 = @(t,x) sqrt(abs(x));
t2 = [0 5];
x2_0 = 0;
[t_1,x_1] = ode45(ode1, t1, x1_0);
[t_2,x_2] = ode45(ode2, t2, x2_0);
figure(1);
plot(t_1,x_1);
 
figure(2);
plot(t_2,x_2);