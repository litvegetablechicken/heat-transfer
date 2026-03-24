function F = double_pipe_residual_L_SA_A_q(x, param)
Nz  = param.geom.Nz;
NzL = param.geom.NzL;
NzSA = param.geom.NzSA;
NzA = param.geom.NzA;

q_d = x(1:Nz);
q_do = x(Nz + 1 : 2*Nz);
q_d_L = q_d(1:NzL);
q_d_SA = q_d(NzL+1 : NzL + NzSA);
q_d_A  = q_d(NzL + NzSA + 1: NzL + NzSA + NzA);
q_do_L = q_do(1:NzL);
q_do_SA = q_do(NzL+1 : NzL + NzSA);
q_do_A  = q_d(NzL + NzSA + 1: NzL + NzSA + NzA);

L_L  = x(end-2);
L_SA = x(end-1);
L_A = x(end);

%% Liquid
out_liquid = Liquid_q(param, q_d_L, L_L);
out_wall_L   = wall_q(param, q_d_L, q_do_L);
outE_L       = external_tube_q(param, q_do_L, L_L,NzL);

r_bp  = out_liquid.T(end) - param.liquid.Tbp;              % 1x1
r_dT_L  = out_wall_L.dT_diff;                                   % Nz x 1
r_dT2_L = (out_liquid.Tw - outE_L.Tw) - out_wall_L.dT;            % Nz x 1
%% SA
Boundary.SA.G_L = out_liquid.G;
Boundary.SA.G_V = 0;
param.external.T_in_ex = outE_L.T_ex(end); % K
out_SA = sub_annular_q(param,Boundary,q_d_SA,L_SA);
out_wall_SA   = wall_q(param, q_d_SA, q_do_SA);
outE_SA       = external_tube_q(param, q_do_SA, L_SA, NzSA);

r_term_SA = out_SA.term;
r_dT_SA  = out_wall_SA.dT_diff;                                   % Nz x 1
r_dT2_SA = (out_SA.Tw - outE_SA.Tw) - out_wall_SA.dT;            % Nz x 1

%% A
Boundary.A.G_L = out_SA.G_L(end);
Boundary.A.G_V = out_SA.G_V(end);
out_A = annular_heat_transfer(param, q_d_A,Boundary,L_A);
out_wall_A   = wall_q(param, q_d_A, q_do_A);
outE_A       = external_tube_q(param, q_do_A, L_A,NzA);
r_dT_A  = out_wall_A.dT_diff;                                   % Nz x 1
r_dT2_A = (out_A.Tw - outE_A.Tw) - out_wall_A.dT;            % Nz x 1
r_term_A = out_A.term;

%% 总误差
F = [[r_dT_L;r_dT_SA;r_dT_A];[r_dT2_L;r_dT2_SA;r_dT2_A];r_bp;r_term_SA;r_term_A];
end
