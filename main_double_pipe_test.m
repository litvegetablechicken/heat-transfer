clear; clc; close all

%% ===== 0) 参数与网格设置 =====
Di = 7.9E-3;
Do = 9.5E-3;
De = 16.0E-3;


kw = 385.0;
G  = 289.0;
G_ex = 211.0; 
NzL = 101;          % Liquid 段网格数
NzSA = 101*2;          % sub_annular 段网格数
% Nz   = NzL + NzSA;      % 总网格数（壁温/外侧共用）
Nz = NzL;
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
param.liquid.heat_capacity= 1.2749;         % kJ/kgK
param.liquid.T_in         = 22 + 273.15;  % K

% --- external/annulus side ---
param.external.G_ex = G_ex;
param.external.m_ex = param.external.G_ex * param.geom.A_ex;

param.external.density_ex       = 989.976875;
param.external.V                = param.external.m_ex / param.external.density_ex;
param.external.viscosity_ex     = 0.000610821;
param.external.thermal_cond_ex  = 0.633815;
param.external.heat_capacity_ex = 4.18095625;  % kJ/kgK
param.external.dir              = 1;
param.external.energy_balance_cal = "energy_bal";
param.external.T_in_ex          = 55.5+273.15; % K
% --- sub_annular side ---

param.sub_annular.P       = 11.5;   % 系统压力
param.sub_annular.Pc      = 43.96;   % 临界压力
param.sub_annular.Mw      = 86.48;   % 分子量
param.sub_annular.DH_vap  = 178.83*1E3;   % 汽化潜热 J/kg
param.sub_annular.g       = 9.81;
param.sub_annular.rhoL    = 1175.4;   % 液相密度
param.sub_annular.rhoV    = 49.148;   % 汽相密度
%% ===== 1) fsolve 初值 =====
Tw_guess = (param.liquid.T_in + param.external.T_in_ex)/2;
Tw_d0  = 320.467024683289 * ones(Nz,1);
Tw_do0 = 320.501859537161 * ones(Nz,1);

L10 = 0.05;           % Liquid 段初猜长度
L20 = 0.25;           % sub_annular 段初猜长度

x0 = [Tw_d0; Tw_do0; L10; L20];
% x0 = [Tw_d0; Tw_do0; L10];
opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-15, ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',2e5, ...
    'UseParallel',true);

[xsol,~,exitflag,output] = fsolve(@(x) double_pipe_L_SA_test(x,param), x0, opts);

disp(output.message); fprintf("exitflag=%d\n",exitflag);

%% ===== 2) 解包 =====
Tw_d  = xsol(1:Nz);
Tw_do = xsol(Nz+1:2*Nz);
L1    = xsol(end-1);
L2    = xsol(end);
% L1 = xsol(end);
% Ltot  = L1 ;

%% ===== 3) 后处理：三模块重新算一遍 =====
% 内侧分段
Tw_d1  = Tw_d(1:NzL);
Tw_d2  = Tw_d(NzL+1:end);
Tw_d2=[];
outL1 = Liquid(param, Tw_d1, L1);
outL2 = sub_annular(param, Tw_d2, L2);

% 拼成总内侧 q_in 和 T（便于画图）
z = outL1.z;


q_in = [outL1.q];
T_in = [outL1.T];    % sub_annular 是常数Tbp

% 外侧与壁（用总长度Ltot、总网格N）
outE = external_tube(param, Tw_do, Ltot);
outW = wall(param, Tw_d, Tw_do);

%% ===== 4) 画图 =====
figure;
plot(z, T_in-273.15, z, outE.T_ex-273.15, z, Tw_d-273.15, z, Tw_do-273.15);
grid on; xlabel('z (m)'); ylabel('T (°C)');
legend('T_{in}','T_{ex}','T_{w,d}','T_{w,do}');

figure;
plot(z, q_in, z, outW.q_d); grid on
xlabel('z (m)'); ylabel('q'); legend('q_{in}','q_d');

figure;
plot(z, outE.q_ex, z, outW.q_do); grid on
xlabel('z (m)'); ylabel('q'); legend('q_{ex}','q_{do}');
