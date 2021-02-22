function dydt = odefcn(t,y,V)
%H -> 1
%I -> 2
%Z -> 3
%D -> 4
dydt = zeros(4,1);
dydt(1) = V.b*y(1) - V.b_d*y(1)*y(1) - V.d*y(1) - V.i*y(1)*y(3);
dydt(2) = -V.a*y(2) + V.i*y(1)*y(3) - V.d*y(2);
dydt(3) = V.a*y(2) + V.r * y(4) - V.n*y(3)*y(1);
dydt(4) = V.d*y(1) + V.d*y(2) - V.r*y(4) + V.n*y(3)*y(1);