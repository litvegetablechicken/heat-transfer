function F = double_pipe(x, param)
Nz1 = param.geom.NzL;
Nz2 = param.geom.NzSA;
N   = param.geom.Nz;
Tw_d  = x(1:N);
Tw_do = x(N+1:2*N);
L_liquid    = x(end-1);
L_SA    = x(end);
Ltot  = L_liquid + L_SA;
% ---- 分段壁温 ----
Tw_d_liquid = Tw_d(1:Nz1);
Tw_d_SA = Tw_d(Nz1+1:end);
% ---- wall(总网格) ----
outW = wall(param, Tw_d, Tw_do);
% ---- inner: Liquid + sub_annular ----
out_liquid = Liquid(param, Tw_d_liquid, L_liquid);
out_SA = sub_annular(param, Tw_d_SA, L_SA);
q_in = [out_liquid.q; out_SA.q];
outE = external_tube(param, Tw_do, Ltot);
% ---- 残差 ----
r_inner = q_in      - outW.q_d;     % N×1
r_outer = outE.q_ex - outW.q_do;    % N×1
% 段间条件：Liquid 段末端到 Tbp
r_bp   = out_liquid.T(end) - param.liquid.Tbp;

% 终止条件：sub_annular 末端 u_ast_C = 1
r_term = out_SA.term;  % outL2.term = u_ast_C(end)-1
F = [r_inner; r_outer; r_bp; r_term];
end