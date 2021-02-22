function dydt = odefcn(t,y,V)
%H -> 1
%I -> 2
%Z -> 3
%D -> 4
%Q -> 5
dydt = zeros(5,1);
dydt(1) = V.b*y(1) - V.b_d*y(1)*y(1) - V.d*y(1) - V.i*y(1)*y(3);
dydt(2) = -V.a*y(2) + V.i*y(1)*y(3) - V.d*y(2) - V.q_i*y(2);
dydt(3) = V.a*y(2) + V.r * y(4) - V.n*y(3)*y(1) - V.q_z*y(3);
dydt(4) = V.d*y(1) + V.d*y(2) - V.r*y(4) + V.n*y(3)*y(1) + V.d_q*y(5);
dydt(5) = V.q_i*y(2) + V.q_z*y(3) - V.d_q*y(5);
