%IRKTemplate(ButcherArray, f, dfdx, T, x0)
clc;
close all;
clear all;
A2 = [(1/4)  (1/4-(sqrt(3)/6)); 
     (1/4 + (sqrt(3)/6)) (1/4)];
c2 = [(1/2 - (sqrt(3)/6)); 
     (1/2 + (sqrt(3)/6))];
b2 = [(1/2); 
      (1/2)];
ButcherExample2 = struct('A',A2,'b',b2,'c',c2);

lambda = -1000000;
t_final = 10000;
t_step = 0.4;
x0 = 1;

f = @(t,x) lambda*x;
dfdx = @(t,x) lambda; 

dt = 0.4;
T = 0:dt:t_final;

S = IRKTemplate(ButcherExample2, f, dfdx, T, x0);


plot(S)