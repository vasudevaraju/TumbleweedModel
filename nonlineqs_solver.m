function [v,fval] = nonlineqs_solver(ts, lhs, states, rm,FW, g, Ixx, Iyy, Izz, m, mf, mu0, R, d, v_thres, minr, maxr)
u0 = [d/2;d/2;d/2];
%v = fsolve(@(u)nonlineqs(u,ts, lhs, states, rm, g, Ixx, Iyy, Izz, m, mf, mu0, R, d, v_thres),u0);
%[v, fval] = fmincon(@(u)nonlineqs(u,ts, lhs, states, rm, g, Ixx, Iyy, Izz, m, mf, mu0, R, d, v_thres),u0,[-1,0,0;0,-1,0;0,0,-1;1,0,0;0,1,0;0,0,1],[-0.2;-0.2;-0.2;0.45;0.45;0.45]);
[v, fval] = fminunc(@(u)nonlineqs(u,ts, lhs, states, rm, FW,g, Ixx, Iyy, Izz, m, mf, mu0, R, d, v_thres, minr, maxr),u0);
end

function eqs = nonlineqs(u, ts, lhs, states, rm,FW, g, Ixx, Iyy, Izz, m, mf, mu0, R, d, v_thres,minr, maxr)
I_T = [Ixx,0,0;0,Iyy,0;0,0,Izz];
r_m1T = [u(1);0;0];             T_v_m1T = [rm(7);0;0]; T_a_m1T = [rm(13);0;0];
r_m2T = [-(d - u(1));0;0];      T_v_m2T = [rm(8);0;0]; T_a_m2T = [rm(14);0;0];
r_m3T = [0;u(2);0];             T_v_m3T = [0;rm(9);0]; T_a_m3T = [0;rm(15);0];
r_m4T = [0;-(d - u(2));0];      T_v_m4T = [0;rm(10);0]; T_a_m4T = [0;rm(16);0];
r_m5T = [0;0;u(3)];             T_v_m5T = [0;0;rm(11)]; T_a_m5T = [0;0;rm(17)];
r_m6T = [0;0;-(d - u(3))];      T_v_m6T = [0;0;rm(12)]; T_a_m6T = [0;0;rm(18)];
r_CMT = m*(r_m1T + r_m2T + r_m3T + r_m4T + r_m5T + r_m6T)/(6*m+mf);
T_v_CMT = m*(T_v_m1T + T_v_m2T + T_v_m3T + T_v_m4T + T_v_m5T + T_v_m6T)/(6*m+mf);
T_a_CMT = m*(T_a_m1T + T_a_m2T + T_a_m3T + T_a_m4T + T_a_m5T + T_a_m6T)/(6*m+mf);
r_CMTx = Vx(r_CMT);r_m1Tx = Vx(r_m1T);r_m2Tx = Vx(r_m2T);r_m3Tx = Vx(r_m3T);
r_m4Tx = Vx(r_m4T);r_m5Tx = Vx(r_m5T);r_m6Tx = Vx(r_m6T);
q0 = states(10);q1=states(11);q2=states(12);q3=states(13);
T_c_O = [q0^2 + q1^2 - q2^2 - q3^2, 2*(q1*q2 + q0*q3), 2*(q1*q3 - q0*q2);...
       2*(q1*q2 - q0*q3), q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 + q0*q1);...
       2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1), q0^2 - q1^2 - q2^2 + q3^2];
O_c_T = transpose(T_c_O);  
r_PT = [0;0;-R];
r_PTx = Vx([0;0;-R]);
O_w_T = [states(1);states(2);states(3)]; O_w_Tx =Vx(O_w_T);

O_vel_TO = Vx(O_c_T*O_w_T)*[0;0;R];
u_TO = O_vel_TO(1);
v_TO = O_vel_TO(2);

murr = mu0*tanh((u_TO^2 + v_TO^2)/v_thres); 
roll_rest_coeff = [1,0,-murr*u_TO/sqrt(u_TO^2 + v_TO^2);0,1,-murr*v_TO/sqrt(u_TO^2 + v_TO^2);0,0,1];
FG = [0;0;-(mf+6*m)*g];
FR  = -murr*norm(FG)/sqrt(u_TO^2 + v_TO^2)*[u_TO; v_TO; 0];

Amat = I_T + (mf + 6*m)*r_CMTx*T_c_O*r_PTx*O_c_T - m*r_m1Tx*r_m1Tx - m*r_m2Tx*r_m2Tx...
       - m*r_m3Tx*r_m3Tx - m*r_m4Tx*r_m4Tx - m*r_m5Tx*r_m5Tx - m*r_m6Tx*r_m6Tx...
       - (mf + 6*m)*Vx(T_c_O*r_PT)*Vx(T_c_O*r_PT) + (mf + 6*m)*Vx(T_c_O*r_PT)*r_CMTx;
  %some mistakes here. check again 
Bvec = -cross(O_w_T,I_T*O_w_T) - m*cross(r_m1T,2*O_w_Tx*T_v_m1T + O_w_Tx*O_w_Tx*r_m1T)...
       - m*r_m2Tx*(T_a_m2T + 2*O_w_Tx*T_v_m2T + O_w_Tx*O_w_Tx*r_m2T)...
       - m*r_m3Tx*(T_a_m3T + 2*O_w_Tx*T_v_m3T + O_w_Tx*O_w_Tx*r_m3T)...
       - m*r_m4Tx*(T_a_m4T + 2*O_w_Tx*T_v_m4T + O_w_Tx*O_w_Tx*r_m4T)...
       - m*r_m5Tx*(T_a_m5T + 2*O_w_Tx*T_v_m5T + O_w_Tx*O_w_Tx*r_m5T)...
       - m*r_m6Tx*(T_a_m6T + 2*O_w_Tx*T_v_m6T + O_w_Tx*O_w_Tx*r_m6T)...
       + r_CMTx*T_c_O*FG + (mf + 6*m)*cross(T_c_O*r_PT, T_a_CMT + ...
         2*O_w_Tx*T_v_CMT + O_w_Tx*O_w_Tx*r_CMT) - T_c_O*cross(r_PT, FG) ...
         - T_c_O*cross(r_PT, FW);%- T_c_O*cross(r_PT, FR);
eqs = norm(inv(Amat)*Bvec - lhs) + 10*max([0,(u(1) -maxr)])^2 + 10*max([0,(minr - u(1))])^2 + 10*max([0,(u(2) -maxr)])^2 + 10*max([0,(minr - u(2))])^2 + 10*max([0,(u(3) -maxr)])^2 + 10*max([0,(minr - u(3))])^2;
end