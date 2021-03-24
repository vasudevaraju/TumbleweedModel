function [v, fval] = controlsolver(ts, lhs, states, rm,FW, g, Ixx, Iyy, Izz, m, mf, mu0, R, d, v_thres, minr, maxr)

coder.extrinsic('nonlineqs_solver');
[v,  fval] = nonlineqs_solver(ts, lhs, states, rm,FW, g, Ixx, Iyy, Izz, m, mf, mu0, R, d, v_thres, minr, maxr);
end