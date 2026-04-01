function F = double_pipe_residual_vapor(x, param,L,T_in)
NzV = param.geom.NzV;
q_d = x(1:NzV);
q_do = x(NzV + 1 : 2*NzV);


%% Liquid
out_liquid = Vapor_q(param, q_d, L,T_in);
out_wall_V   = wall_q(param, q_d, q_do);
outE_V       = external_tube_q(param, q_do, L,NzV);

r_dT_L  = out_wall_V.dT_diff;                                   % Nz x 1
r_dT2_L = (out_liquid.Tw - outE_V.Tw) - out_wall_V.dT;            % Nz x 1

%% 总误差
F = [r_dT_L;r_dT2_L];
end
