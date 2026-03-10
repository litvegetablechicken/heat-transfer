function F = double_pipe_residual_SA_q(x, param,Boundary)
NzSA = param.geom.NzSA;

q_d  = x(1:NzSA);
q_do = x(NzSA+1 : NzSA+NzSA);
L_SA = x(end);
out_SA = sub_annular_q(param,Boundary,q_d,L_SA);
out_wall   = wall_q(param, q_d, q_do);
outE       = external_tube_q(param, q_do, L_SA);

r_dT_w  = out_wall.dT_diff;                                   % Nz x 1
r_dT2 = (out_SA.Tw - outE.Tw) - out_wall.dT;            % Nz x 1
r_term = out_SA.term;
F = [r_dT_w; r_dT2;r_term];
end
