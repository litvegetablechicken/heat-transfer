function F = double_pipe_mist_q(x, param,Boundary)
NzM = param.geom.NzM;

q_d  = x(1:NzM);
q_do = x(NzM+1 : NzM+NzM);
L_M = x(end);
out_M = mist_q(param, q_d,Boundary,L_M);
out_wall   = wall_q(param, q_d, q_do);
outE       = external_tube_q(param, q_do, L_M,NzM);

r_dT_w  = out_wall.dT_diff;                                   % Nz x 1
r_dT2 = (out_M.Tw - outE.Tw) - out_wall.dT;            % Nz x 1
r_term = out_M.G_ED(end);
F = [r_dT_w; r_dT2;r_term];
end