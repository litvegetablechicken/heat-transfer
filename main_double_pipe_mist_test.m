clear; clc; close all

%% ===== 0) 参数与网格设置 =====
name = 'R22';
% geometery
Di = 7.9E-3;
Do = 9.5E-3;
De = 16.0E-3;
% wall
kw = 385.0;
% fluid
P = 11.5;
Pc = 43.96;
Mw = 86.48;
DH_vap = 178.83*1E3; % J/kg
g = 9.81;
rhoL = 1175.4;
rhoV = 49.148;
muL  = 0.00015749;      % Pa*s
kL   = 0.081729; % W/m/K  (占位)
cpL  = 1.2749*1e3; % J/kg/K
muV  = 0.000012712;          % Pa*s  
kV   = 0.011696;            % W/m/K  (占位)
cpV  = 0.89952*1e3;           % J/kg/K  (占位)
sigma = 0.0075698;     % N/m (占位)

rhoex = 989.976875; % kg/m3
muex = 0.000610821; %  Pas 
kex = 0.633815; % w/m/K
cpex = 4.18095625*1000; % J/kg
A_ant = 9.357495827;     % Antoine A (log10(P_bar)=A - B/(C+T_C))
B_ant = 945.3909755;  % Antoine B
C_ant = -14.964;   % Antoine C
Tbp = 28.6 + 273.15;

% flow
G  = 289.0;
G_ex = 211.0; 
T_in = 22 + 273.15; %K
T_in_ex = 305.39435; % K
dir              = 1; % 1 顺流 -1 逆流

NzM = 1001*2;
Nz = NzM;
%% ========= 0) 统一参数 param =========
param = struct();

param.geom.NzM = NzM;


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
param.liquid.Tbp = Tbp;
param.liquid.m            = param.liquid.G * param.geom.A_i;
param.liquid.density      = rhoL;
param.liquid.viscosity    = muL;
param.liquid.thermal_cond = kL;
param.liquid.heat_capacity= cpL;         % J/kgK
param.liquid.T_in         = T_in;  % K

% --- external/annulus side ---
param.external.G_ex = G_ex;
param.external.m_ex = param.external.G_ex * param.geom.A_ex;

param.external.density_ex       = rhoex;
param.external.V                = param.external.m_ex / param.external.density_ex;
param.external.viscosity_ex     = muex;
param.external.thermal_cond_ex  = kex;
param.external.heat_capacity_ex = cpex;  % kJ/kgK
param.external.dir              = dir; % 1 顺流 -1 逆流
param.external.energy_balance_cal = "energy_bal";
param.external.T_in_ex          = T_in_ex; % K
% --- sub_annular side ---

param.sub_annular.P       = P;   % 系统压力
param.sub_annular.Pc      = Pc;   % 临界压力
param.sub_annular.Mw      = Mw;   % 分子量
param.sub_annular.DH_vap  = DH_vap;   % 汽化潜热 J/kg
param.sub_annular.g       = g;
param.sub_annular.rhoL    = rhoL;   % 液相密度
param.sub_annular.rhoV    = rhoV;   % 汽相密度

% annular
%  2  Fluid properties (example: R245fa)
param.model.delta_calculation        = true;
param.model.h_calculation            = "correlation_0";
param.model.calculate_entraintment_0 = true;

param.fluid.name = name;
param.fluid.Tbp = Tbp;
param.fluid.rhoL = rhoL;          % kg/m3
param.fluid.muL  = muL;      % Pa*s
param.fluid.kL   = kL;        % W/m/K
param.fluid.cpL  = cpL;        % J/kg/K ✅

param.fluid.rhoV = rhoV;          % kg/m3
param.fluid.muV  = muV;          % Pa*s  (占位，建议换真实)
param.fluid.kV   = kV;            % W/m/K  (占位)
param.fluid.cpV  = cpV;           % J/kg/K  (占位)

param.fluid.rhoED = param.fluid.rhoL;
param.fluid.sigma = sigma;     % N/m (占位)
param.fluid.DH_vap = DH_vap;          % J/kg ✅
param.fluid.g = g;

% Cooper 参数（用 Pa，P/Pc 无量纲）
param.fluid.P  = P;    % Pa
param.fluid.Pc = Pc;   % Pa
param.fluid.Mw = 86.48;     % 分子量（用于幂次）
param.fluid.A_ant = A_ant;     % Antoine A (log10(P_bar)=A - B/(C+T_C))
param.fluid.B_ant = B_ant;  % Antoine B
param.fluid.C_ant = C_ant;   % Antoine C

Boundary.Mist.G_L0 = 0;
Boundary.Mist.G_V0 = 291.0839;
Boundary.Mist.G_ED0 = 0.5005622;

x0 = init_double_pipe_q(param,T_in_ex,Tbp); % 1-Nz 是q_d，Nz+1-2Nz是q_do，2Nz+1 是L_L, 2Nz+2是L_SA 


opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-20, ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',2e5,...
    UseParallel=true);

[xsol,feval,exitflag,output] = fsolve(@(x)double_pipe_mist_q(x, param,Boundary), x0, opts);

q_d  = xsol(1:NzM);
q_do = xsol(NzM+1 : NzM+NzM);
L_M = xsol(end);
out_M = mist_q(param, q_d,Boundary,L_M);
out_wall   = wall_q(param, q_d, q_do);
outE       = external_tube_q(param, q_do, L_M,NzM);