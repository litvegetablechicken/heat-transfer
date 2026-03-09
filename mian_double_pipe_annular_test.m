clear; clc; close all
%TO DO anular模型中换热系数作为自变量迭代
%% ===== 0) 参数与网格设置 =====
Di = 7.9E-3;
Do = 9.5E-3;
De = 16.0E-3;

kw   = 385.0;
G_ex = 211.0;

Nz = 2;   % annular 段网格数（单独测试 annular）

param = struct();

% --- geom ---
param.geom.Nz = Nz;
param.geom.Di = Di;
param.geom.Do = Do;
param.geom.De = De;
param.geom.A_i  = pi*Di^2/4;
param.geom.A_ex = pi*(De^2 - Do^2)/4;
param.geom.Dh   = De - Do;

% --- wall ---
param.wall.kw = kw;

% --- external side（全部 SI） ---
param.external.G_ex = G_ex;
param.external.m_ex = param.external.G_ex * param.geom.A_ex;
param.external.density_ex       = 989.976875;
param.external.V                = param.external.m_ex / param.external.density_ex;
param.external.viscosity_ex     = 0.000610821;
param.external.thermal_cond_ex  = 0.633815;
param.external.heat_capacity_ex = 4.18095625e3;   % J/kg/K  ✅
param.external.dir              = 1;
param.external.energy_balance_cal = "energy_bal";
param.external.T_in_ex          = 322.6129315;

% --- saturation temperature ---
param.liquid.Tbp = 28.6 + 273.15;

% --- annular switches ---
param.annular.delta_calculation        = "triangular_rel";
param.annular.entraintment_calculation = "cal_entraintment";
param.annular.h_calculation            = "correlation_0";

% --- annular properties (SI) ---
param.annular.rhoL = 1175.4;          % kg/m3
param.annular.muL  = 0.00015749;      % Pa*s
param.annular.kL   = 0.081729;        % W/m/K
param.annular.cpL  = 1.2749*1e3;        % J/kg/K ✅

param.annular.rhoV = 49.148;          % kg/m3
param.annular.muV  = 0.000012712;          % Pa*s  (占位，建议换真实)
param.annular.kV   = 0.011696;            % W/m/K  (占位)
param.annular.cpV  = 0.89952*1e3;           % J/kg/K  (占位)

param.annular.rhoED = param.annular.rhoL;
param.annular.surface_tension = 0.0075698;     % N/m (占位)
param.annular.DH_vap = 178.83e3;          % J/kg ✅
param.annular.g = 9.81;

% Cooper 参数（用 Pa，P/Pc 无量纲）
param.annular.P  = 11.5e5;    % Pa
param.annular.Pc = 43.96e5;   % Pa
param.annular.Mw = 86.48;     % 分子量（用于幂次）
param.annular.A_ant = 4.0;     % Antoine A (log10(P_bar)=A - B/(C+T_C))
param.annular.B_ant = 1000.0;  % Antoine B
param.annular.C_ant = 200.0;   % Antoine C
% inlet mass flux (kg/m2/s)
param.annular.G_L0 = 223.503399116298;
param.annular.G_V0 = 65.4966008837023;

state0 = struct();
state0.G_L0 = param.annular.G_L0;
state0.G_V0 = param.annular.G_V0;
state0.calculate_entraintment_0 = true;

%% ===== 1) fsolve 初值 =====
Tw_d0  = (303.5629929) * ones(Nz,1);
Tw_do0 = (303.644088370449) * ones(Nz,1);
% L0 = 1;
param.L = 0.00166667;                     % m
delta0 = 0.00060943 * ones(Nz,1);     % ✅ 数组初猜 (m)
% 物理范围建议：0 < delta < Di/2
% 你也可以给个线性分布作为初猜

% x0 = [Tw_d0; Tw_do0; L0; delta0];
x0 = [Tw_d0; Tw_do0; delta0];

opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-20, ...
    'StepTolerance',1e-20, ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',2e5,...
    UseParallel=false);

[xsol,feval,exitflag,output] = fsolve(@(x) double_pipe_residual_annular(x, param, state0), x0, opts);

disp(output.message); fprintf("exitflag=%d\n",exitflag);

%% ===== 2) 解包并后处理 =====
N = Nz;

Tw_d  = xsol(1:N);
Tw_do = xsol(N+1:2*N);
% L_ann = xsol(2*N+1);
% delta = xsol(2*N+2:end);
delta = xsol(2*N+1:end);
outA = annular(param, Tw_d, param.L, delta, state0);
outE = external_tube(param, Tw_do, param.L);
outW = wall(param, Tw_d, Tw_do);

z = outA.z;

%% ===== 3) 画图 =====


figure;
plot(z, outA.G_L, z, outA.G_V, z, outA.G_ED); grid on
xlabel('z (m)'); ylabel('G (kg/m^2/s)'); legend('G_L','G_V','G_{ED}');

figure;
plot(z, Tw_d, z, Tw_do); grid on
xlabel('z (m)'); ylabel('T (K)'); legend('Tw_d','Tw_do');
title('delta consistency: in vs model');
