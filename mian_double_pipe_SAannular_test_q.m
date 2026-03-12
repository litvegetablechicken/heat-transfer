clear; clc; close all

%% ===== 0) 参数与网格设置 =====
Di = 7.9E-3;
Do = 9.5E-3;
De = 16.0E-3;
P = 11.5;
Pc = 43.96;
Mw = 86.48;
DH_vap = 178.83*1E3;
g = 9.81;
rhoL = 1175.4;
rhoV = 49.148;
kw = 385.0;
G  = 289.0;
G_ex = 211.0; 
NzL = 1001;          % Liquid 段网格数
NzSA = 1001;          % sub_annular 段网格数
% Nz   = NzL + NzSA;      % 总网格数（壁温/外侧共用）
Nz = NzL + NzSA;
%% ========= 0) 统一参数 param =========
param = struct();

param.geom.NzL = NzL;
param.geom.NzSA = NzSA;
param.geom.Nz   = Nz;

param.geom.Di = Di;
param.geom.Do = Do;
param.geom.De = De;

param.geom.A_i  = pi*param.geom.Di^2/4;
param.geom.A_ex = pi*(param.geom.De^2 - param.geom.Do^2)/4;
param.geom.Dh   = param.geom.De - param.geom.Do;
% --- wall ---
param.wall.kw = kw;  % thermal_cond_wall

% --- liquid side ---
param.liquid.G  = G; % 你原来的 G（用于算 m = G*A_i）
param.liquid.Tbp = 28.6+273.15;
param.liquid.m            = param.liquid.G * param.geom.A_i;
param.liquid.density      = 1175.4;
param.liquid.viscosity    = 0.00015749;
param.liquid.thermal_cond = 0.081729;
param.liquid.heat_capacity= 1.2749*1000;         % J/kgK
param.liquid.T_in         = 22 + 273.15;  % K

% --- external/annulus side ---
param.external.G_ex = G_ex;
param.external.m_ex = param.external.G_ex * param.geom.A_ex;

param.external.density_ex       = 989.976875;
param.external.V                = param.external.m_ex / param.external.density_ex;
param.external.viscosity_ex     = 0.000610821;
param.external.thermal_cond_ex  = 0.633815;
param.external.heat_capacity_ex = 4.18095625*1000;  % kJ/kgK
param.external.dir              = 1; % 1 顺流 -1 逆流
param.external.energy_balance_cal = "energy_bal";
param.external.T_in_ex          = 55.5+273.15; % K
% --- sub_annular side ---

param.sub_annular.P       = P;   % 系统压力
param.sub_annular.Pc      = Pc;   % 临界压力
param.sub_annular.Mw      = Mw;   % 分子量
param.sub_annular.DH_vap  = DH_vap;   % 汽化潜热 J/kg
param.sub_annular.g       = g;
param.sub_annular.rhoL    = rhoL;   % 液相密度
param.sub_annular.rhoV    = rhoV;   % 汽相密度

%% ===== 1) fsolve 初值 =====
x0 = init_double_pipe_q(param); % 1-Nz 是q_d，Nz+1-2Nz是q_do，2Nz+1 是L_L, 2Nz+2是L_SA 

opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-15, ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',2e5, ...
    'UseParallel',true);

[x_sol,f,exitflag,output] = fsolve(@(x) double_pipe_residual_L_SA_q(x, param), x0, opts);

disp(output.message); fprintf("exitflag=%d\n",exitflag);

% 解析结果
q_d = x_sol(1:Nz);
q_do = x_sol(Nz + 1 : 2*Nz);
L_L  = x_sol(end-1);
L_SA = x_sol(end);
L_total = L_L + L_SA;
q_d_L = q_d(1:NzL);
q_d_SA = q_d(NzL+1 : end);
q_do_L = q_do(1:NzL);
q_do_SA = q_do(NzL+1 : end);

%% 后处理
% Liquid
out_liquid = Liquid_q(param, q_d_L, L_L);
out_wall_L   = wall_q(param, q_d_L, q_do_L);
outE_L       = external_tube_q(param, q_do_L, L_L,NzL);
% SA
Boundary.SA.G_L = out_liquid.G;
Boundary.SA.G_V = 0;
out_SA = sub_annular_q(param,Boundary,q_d_SA,L_SA);
out_wall_SA   = wall_q(param, q_d_SA, q_do_SA);
outE_SA       = external_tube_q(param, q_do_SA, L_SA, NzSA);

% z = outE.z;
% figure;
% plot(z,outE.T_ex, z, outE.Tw, [out_liquid.z;out_SA.z + L_L], [out_liquid.Tw;out_SA.Tw], [out_liquid.z;out_SA.z + L_L], [out_liquid.T; out_SA.T])
% legend('Tex','Two','Tw','Torc');
% z = out_SA.z;
% figure;
% plot(z, out_SA.Tw, z, outE.Tw);
% grid on; xlabel('z (m)'); ylabel('T (K)');
% legend('Tw','Two');
% figure;
% plot(z, out_SA.T, z, outE.T_ex);
% grid on; xlabel('z (m)'); ylabel('T (K)');
% legend('Torc','Tex');
