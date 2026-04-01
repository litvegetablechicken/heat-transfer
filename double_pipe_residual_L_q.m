function F = double_pipe_residual_L_q(x, param)

NzL = param.geom.NzL;
q_d = x(1:NzL);
q_do = x(NzL + 1 : 2*NzL);
L_L  = x(end);


%% Liquid
out_liquid = Liquid_q(param, q_d, L_L);
out_wall_L   = wall_q(param, q_d, q_do);
outE_L       = external_tube_q(param, q_do, L_L,NzL);

r_bp  = out_liquid.T(end) - param.liquid.Tbp;              % 1x1
r_dT_L  = out_wall_L.dT_diff;                                   % Nz x 1
r_dT2_L = (out_liquid.Tw - outE_L.Tw) - out_wall_L.dT;            % Nz x 1

%% 总误差
F = [r_dT_L;r_dT2_L;r_bp*10];
end
