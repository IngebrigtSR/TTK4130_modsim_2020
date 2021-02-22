%%Task 1a) and 1b)
syms x y;
f_sym1 = x*y - 2;
f_sym2 = ((x^4)/4) + ((y^3)/3) - 1;
f_sym = [f_sym1;
         f_sym2];
    
F = matlabFunction(f_sym, "vars", {[x;y]});

J_sym = jacobian(f_sym, [x;y]);
J = matlabFunction(J_sym, "vars", {[x;y]});

x0 = [-1;
      -1];
tol = 0.1;
N = 100;

[S_b, infNorm_b] = NewtonsMethodTemplate(F, J, x0);

%plots 1a) 1b)
% figure(1)
% loglog(infNorm_b);
% figure(2)
% semilogy(infNorm_b);
% title('fig 1: infNorm of residuals')


%%Task 1c)
f_c = (x - 1)*(x - 2)*(x - 3) + 1;
F_c = matlabFunction(f_c, "vars", x);
j_c = jacobian(f_c, x);
J_c = matlabFunction(j_c, "vars", x);
x_c_0 = 3;

[S_c, infNorm_c] = NewtonsMethodTemplate(F_c, J_c, x_c_0);

%plots 1c)
% figure(3)
% plot(S_c');
% title('fig 2: Plot of obtained results in 1c)')

%%Task 1d)
syms x_1 x_2;
f_d1 = x_1 - 1 + (cos(x_2)*x_1 + 1)*cos(x_2);
f_d2 = -x_1*sin(x_2)*(cos(x_2)* x_1 + 1);
f_d = [f_d1;
       f_d2];

F_d = matlabFunction(f_d, "vars", {[x_1;x_2]});
j_d = jacobian(f_d, [x_1;x_2]);
J_d = matlabFunction(j_d, "vars", {[x_1;x_2]});

x_d_0 = [1;
         3];

[S_d, infNorm_d] = NewtonsMethodTemplate(F_d, J_d, x_d_0);

%plots 1d)
% figure(6)
% semilogy(infNorm_d);
% title('fig 3: infNorm of residuals 1d)')

%%Task 1e)
f_e = 100*(x_2 - x_1)^2 + (x_1 - 1)^4;
df_e = [diff(f_e, x_1);
        diff(f_e, x_2)];

F_e = matlabFunction(df_e, "vars", {[x_1;x_2]});

j_e = jacobian(df_e, [x_1;x_2]);
J_e = matlabFunction(j_e, "vars", {[x_1;x_2]});

x_e_0 = [10;
         10];

[S_e, infNorm_e] = NewtonsMethodTemplate(F_e, J_e, x_e_0);

%plots 1e)
figure(7)
plot(infNorm_e);
title('')

%kasper_matte1
Y(x) = 7*x^3 + 9 + sinh(4*x);
Y_d = 21 * x^2 + 4*cosh(4*x);
Y_newt = NewtonsMethodTemplate(Y(x), Y_d, -1, 0.1, 10);



