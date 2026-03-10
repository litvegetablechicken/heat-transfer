function F = double_pipe_residual_liquid_q(x, param)

NzL = param.geom.NzL;
Nz  = param.geom.Nz;

% 为了让 fsolve 方程数 = 未知数数，最好要求 NzL == Nz
if NzL ~= Nz
    error("double_pipe_residual_q:GridMismatch", ...
        "For this version, require param.geom.NzL == param.geom.Nz.");
end

% -----------------------------
% 解析未知量
% x = [q_d(1:NzL); q_do(1:Nz); L]
% -----------------------------
q_d  = x(1:NzL);
q_do = x(NzL+1 : NzL+Nz);
L    = x(end);


% -----------------------------
% 调用三个子模型
% -----------------------------
out_liquid = Liquid_q(param, q_d, L);
out_wall   = wall_q(param, q_d, q_do);
outE       = external_tube_q(param, q_do, L);

% -----------------------------
% 残差
% -----------------------------
r_bp  = out_liquid.T(end) - param.liquid.Tbp;              % 1x1
r_dT  = out_wall.dT_diff;                                   % Nz x 1
r_dT2 = (out_liquid.Tw - outE.Tw) - out_wall.dT;            % Nz x 1

% 总残差
F = [r_bp; r_dT; r_dT2];

end