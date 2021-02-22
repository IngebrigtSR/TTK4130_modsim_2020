%Arbitrary rotational matrix
R = [0.788571   0.377143    0.485714;
 -0.337143   0.925714   -0.171429;
 -0.514286  -0.0285714   0.857143];

%Math
r_00 = trace(R);
T = r_00;
r_11 = R(1,1);
r_22 = R(2,2);
r_33 = R(3,3);
r_vec = [r_00;
        r_11;
        r_22;
        r_33];

syms z_0 z_1 z_2 z_3
z_0 = sqrt(1+2*r_00-T);
z_1 = sqrt(1+2*r_11-T);
z_2 = sqrt(1+2*r_22-T);
z_3 = sqrt(1+2*r_33-T);
z_num = [z_0; z_1; z_2; z_3];
r_ii = max(r_vec);
z_i = abs(sqrt(1+2*r_ii)-r_00);

n_withBigDick = z_0 / 2;
epsilon_i = 0.5 * [z_1; 
                    z_2; 
                    z_3];

%Result
eulerParameters = [n_withBigDick;
                    epsilon_i];


