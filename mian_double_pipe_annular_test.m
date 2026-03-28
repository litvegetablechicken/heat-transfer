clear; clc; close all
%TO DO anular模型中换热系数作为自变量迭代
%% ===== 0) 参数与网格设置 =====
Di = 7.9E-3;
Do = 9.5E-3;
De = 16.0E-3;

kw   = 385.0;
G_ex = 211.0;

NzA = 1001;   % annular 段网格数（单独测试 annular）

param = struct();

% --- geom ---


param.geom.NzA = NzA;
param.geom.Nz  = NzA;
param.geom.Di = Di;
param.geom.Do = Do;
param.geom.De = De;
param.geom.A_i  = pi*Di^2/4;
param.geom.A_ex = pi*(De^2 - Do^2)/4;
param.geom.Dh   = De - Do;

% --- wall ---
param.wall.kw = kw;
% --- liquid side ---
param.liquid.G  = 290; % 你原来的 G（用于算 m = G*A_i）
param.liquid.Tbp = 28.6+273.15;
param.liquid.m            = param.liquid.G * param.geom.A_i;
param.liquid.density     = 1175.4;
param.liquid.viscosity    = 0.00015749;
param.liquid.thermal_cond = 0.081729;
param.liquid.heat_capacity= 1.2749*1000;         % J/kgK
param.liquid.T_in         = 22 + 273.15;  % K
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
param.external.T_in_ex          = 322.61295; % To be modified

% Model options



% --- annular switches ---
param.model.delta_calculation        = true;
param.model.h_calculation            = "correlation_0";
param.model.calculate_entraintment_0 = true;

%  2  Fluid properties (example: R245fa)
param.fluid.name = 'R22';
param.fluid.Tbp = 28.6 + 273.15;
param.fluid.rhoL = 1175.4;          % kg/m3
param.fluid.muL  = 0.00015749;      % Pa*s
param.fluid.kL   = 0.081729;        % W/m/K
param.fluid.cpL  = 1.2749*1e3;        % J/kg/K ✅

param.fluid.rhoV = 49.148;          % kg/m3
param.fluid.muV  = 0.000012712;          % Pa*s  (占位，建议换真实)
param.fluid.kV   = 0.011696;            % W/m/K  (占位)
param.fluid.cpV  = 0.89952*1e3;           % J/kg/K  (占位)

param.fluid.rhoED = param.fluid.rhoL;
param.fluid.sigma = 0.0075698;     % N/m (占位)
param.fluid.DH_vap = 178.83e3;          % J/kg ✅
param.fluid.g = 9.81;

% Cooper 参数（用 Pa，P/Pc 无量纲）
param.fluid.P  = 11.5e5;    % Pa
param.fluid.Pc = 43.96e5;   % Pa
param.fluid.Mw = 86.48;     % 分子量（用于幂次）
param.fluid.A_ant = 4.0;     % Antoine A (log10(P_bar)=A - B/(C+T_C))
param.fluid.B_ant = 1000.0;  % Antoine B
param.fluid.C_ant = 200.0;   % Antoine C
% 5  Boundary conditions

Boundary.A.G_L = 223.5034;

Boundary.A.G_V = 65.4966;


%% ===== 1) fsolve 初值 =====
x0 = init_double_pipe_q(param,param.external.T_in_ex ,param.fluid.Tbp); % 1-Nz 是q_d，Nz+1-2Nz是q_do，2Nz+1 是L_L, 2Nz+2是L_SA 


opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-20, ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',2e5,...
    UseParallel=true);

[xsol,feval,exitflag,output] = fsolve(@(x) double_pipe_residual_A_q(x, param,Boundary), x0, opts);

disp(output.message); fprintf("exitflag=%d\n",exitflag);
q_d = xsol(1:NzA);
q_do = xsol(NzA + 1 : 2*NzA);
L_A  = xsol(end);
out_A = annular_heat_transfer(param, q_d,Boundary,L_A);
out_wall   = wall_q(param, q_d, q_do);
outE       = external_tube_q(param, q_do, L_A,NzA);
%% mist
NzM = 1001;
param.geom.NzM = NzM;
param.geom = rmfield(param.geom, 'NzA');
param.geom.Nz = NzM;
Boundary.Mist.G_L0 = out_A.G_L(end);
Boundary.Mist.G_V0 = out_A.G_V(end);
Boundary.Mist.G_ED0 = out_A.G_ED(end);
param.external.T_in_ex          = outE.T_ex(end); % To be modified

x0 = init_double_pipe_q(param,outE.T_ex(end),param.liquid.Tbp); % 1-Nz 是q_d，Nz+1-2Nz是q_do，2Nz+1 是L_L, 2Nz+2是L_SA 


opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-20, ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',2e5,...
    UseParallel=true);

[xsol,feval,exitflag,output] = fsolve(@(x)double_pipe_mist_q(x, param,Boundary), x0, opts);

q_d_A  = xsol(1:NzM);
q_do_A = xsol(NzM+1 : NzM+NzM);
L_M = xsol(end);
out_M = mist_q(param, q_d_A,Boundary,L_M);
out_wall_M   = wall_q(param, q_d_A, q_do_A);
outE_M      = external_tube_q(param, q_do_A, L_M,NzM);