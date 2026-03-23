function F = double_pipe_residual_L_SA_q(x, param)
Nz  = param.geom.Nz;
NzL = param.geom.NzL;
NzSA = param.geom.NzSA;
q_d = x(1:Nz);
q_do = x(Nz + 1 : 2*Nz);
q_d_L = q_d(1:NzL);
q_d_SA = q_d(NzL+1 : end);
q_do_L = q_do(1:NzL);
q_do_SA = q_do(NzL+1 : end);
L_L  = x(end-1);
L_SA = x(end);
L_total = L_L + L_SA;

%% Liquid
out_liquid = Liquid_q(param, q_d_L, L_L);
out_wall_L   = wall_q(param, q_d_L, q_do_L);
outE_L       = external_tube_q(param, q_do_L, L_L,NzL);

r_bp  = out_liquid.T(end) - param.fluid.Tbp;              % 1x1
r_dT_L  = out_wall_L.dT_diff;                                   % Nz x 1
r_dT2_L = (out_liquid.Tw - outE_L.Tw) - out_wall_L.dT;            % Nz x 1
%% SA
Boundary.SA.G_L = out_liquid.G;
Boundary.SA.G_V = 0;
param.external.T_in_ex = outE_L.T_ex(end); % K
out_SA = sub_annular_q(param,Boundary,q_d_SA,L_SA);
out_wall_SA   = wall_q(param, q_d_SA, q_do_SA);
outE_SA       = external_tube_q(param, q_do_SA, L_SA, NzSA);

r_term = out_SA.term;
r_dT_SA  = out_wall_SA.dT_diff;                                   % Nz x 1
r_dT2_SA = (out_SA.Tw - outE_SA.Tw) - out_wall_SA.dT;            % Nz x 1

%% 总误差
F = [[r_dT_L;r_dT_SA];[r_dT2_L;r_dT2_SA];r_bp;r_term];
end
