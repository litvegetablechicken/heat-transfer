function F = double_pipe_residual_A_q_for_ini(x, param,Boundary,L_A)
NzA = param.geom.NzA;

q_d  = x(1:NzA);
q_do = x(NzA+1 : NzA+NzA);
delta0 = x(NzA+NzA + 1 : end);
out_A = annular_heat_transfer(param, q_d,Boundary,L_A,delta0);
out_wall   = wall_q(param, q_d, q_do);
outE       = external_tube_q(param, q_do, L_A,NzA);

r_dT_w  = out_wall.dT_diff;                                   % Nz x 1
r_dT2 = (out_A.Tw - outE.Tw) - out_wall.dT;            % Nz x 1
F = [r_dT_w; r_dT2;out_A.dif_delta];
end
