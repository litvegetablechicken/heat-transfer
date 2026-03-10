clear; clc; close all

%% ===== 0) 参数与网格设置 =====
Di = 7.9E-3;
Do = 9.5E-3;
De = 16.0E-3;


kw = 385.0;
G  = 289.0;
G_ex = 211.0; 
NzL = 0;          % Liquid 段网格数
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
param.liquid.heat_capacity= 1.2749e3;         % J/kgK
param.liquid.T_in         = 295.15;  % K

% --- external/annulus side ---
param.external.G_ex = G_ex;
param.external.m_ex = param.external.G_ex * param.geom.A_ex;

param.external.density_ex       = 989.976875;
param.external.V                = param.external.m_ex / param.external.density_ex;
param.external.viscosity_ex     = 0.000610821;
param.external.thermal_cond_ex  = 0.633815;
param.external.heat_capacity_ex = 4.18095625e3;  % J/kgK
param.external.dir              = 1;
param.external.energy_balance_cal = "energy_bal";
param.external.T_in_ex          = 327.6121; % K
% --- sub_annular side ---

param.sub_annular.P       = 11.5;   % 系统压力
param.sub_annular.Pc      = 43.96;   % 临界压力
param.sub_annular.Mw      = 86.48;   % 分子量
param.sub_annular.DH_vap  = 178.83*1E3;   % 汽化潜热 J/kg
param.sub_annular.g       = 9.81;
param.sub_annular.rhoL    = 1175.4;   % 液相密度
param.sub_annular.rhoV    = 49.148;   % 汽相密度
Boundary.SA.G_L = 289.0;
Boundary.SA.G_V = 0;
%% ===== 1) fsolve 初值 =====
x0 = init_double_pipe_q(param);

opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-15, ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',2e5, ...
    'UseParallel',true);

[x_sol,~,exitflag,output] = fsolve(@(x) double_pipe_residual_SA_q(x, param,Boundary), x0, opts);

disp(output.message); fprintf("exitflag=%d\n",exitflag);

% 解析结果
q_d_sol  = x_sol(1:Nz);
q_do_sol = x_sol(Nz+1:2*Nz);
L_sol    = x_sol(end);

% 后处理
out_SA = sub_annular_q(param,Boundary,q_d_sol,L_sol);
out_wall   = wall_q(param, q_d_sol, q_do_sol);
outE       = external_tube_q(param, q_do_sol, L_sol);
z = out_SA.z;
figure;
plot(z, out_SA.Tw, z, outE.Tw);
grid on; xlabel('z (m)'); ylabel('T (K)');
legend('Tw','Two');
figure;
plot(z, out_SA.T, z, outE.T_ex);
grid on; xlabel('z (m)'); ylabel('T (K)');
legend('Torc','Tex');
