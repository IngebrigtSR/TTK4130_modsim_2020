syms H(t) I(t) Z(t) D(t) Q(t);

a=1.4*10^(-6);
b=3.1*10^(-8);
b_d=5.6*10^(-16);
d=2.8*10^(-8);
i=2.6*10^(-6);
n=1.4*10^(-6);
r=2.8*10^(-7);
q_i=2.7*10^(-6);
q_z=2.7*10^(-6);
d_q=2.8*10^(-5);
dH_s = H(0)==((b-d)/b_d);
dI_s = I(0)==0;
dZ_s = Z(0)==0;
dD_s = D(0)==0;
dQ_s = Q(0)==0;

timespan = [0 1000];

dH = diff(H) == b*H - b_d*H^2 -d*H -i*H*Z;
dI = diff(I) ==i*H*Z -a*I -d*I-I*q_i;
dZ= diff(Z) == a*I + r*D -n*Z*H-Z*q_z;
dD = diff(D) == d*H + d*I -r*D + n*Z*H+d_q*Q;
dQ = diff(Q) == I*q_i + Z*q_z - d_q*Q;
dAll = diff(H)+diff(I)+diff(Z)+diff(D)+diff(Q) == b*H - b_d*H^2;

odes = [dH; dI; dZ; dD; dQ];
conds = [dH_s;dI_s;dZ_s;dD_s;dQ_s];
[timespan, H(t)] = ode45(diff(H), timespan, dH_s);
[timespan, I(t)] = ode45(diff(I), timespan, dI_s);
[timespan, Z(t)] = ode45(diff(Z), timespan, dZ_s);
[timespan, D(t)] = ode45(diff(D), timespan, dD_s);
[timespan, Q(t)] = ode45(diff(Q), timespan, dQ_s);
plot(timespan, H);

%[Hsol(t),Isol(t),Zsol(t),Dsol(t),Dsol(t),Qsol(t)] =