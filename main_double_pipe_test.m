clear; clc; close all
%% ========= 0) 设定参数 =========
Nz = 101;
% L  = 0.2961;

Di = 7.9E-3;
Do = 9.5E-3;
De = 16.0E-3;


kw = 385.0;
G  = 289.0;
G_ex = 211.0;
%% ========= 0) 统一参数 param =========
param = struct();

% --- 网格/几何 ---
param.geom.NzL = Nz;


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
param.liquid.Tbp = 22+273.15;
param.liquid.m            = param.liquid.G * param.geom.A_i;
param.liquid.density      = 1175.4;
param.liquid.viscosity    = 0.00015749;
param.liquid.thermal_cond = 0.081729;
param.liquid.heat_capacity= 1.2749;         % kJ/kgK
param.liquid.T_in         = 28.6 + 273.15;  % K

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

%% ========= 1) fsolve 初值 =========
Nz = param.geom.NzL;
param.geom.Nz = param.geom.NzL;
Tw_d0  = linspace((param.liquid.T_in + param.external.T_in_ex)/2, ...
                  (param.liquid.T_in + param.external.T_in_ex)-1, Nz).';

Tw_do0 = linspace((param.liquid.T_in + param.external.T_in_ex)/2, ...
                  (param.liquid.T_in + param.external.T_in_ex)-1, Nz).';
L0 = 1;
x0 = [Tw_d0; Tw_do0;L0];

%% ========= 2) fsolve =========
opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-15, ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',2e5,'UseParallel',true);

[xsol, fval, exitflag, output] = fsolve(@(x) double_pipe_residual_liquid(x, param), x0, opts);

fprintf("\n==== fsolve finished ====\n");
fprintf("exitflag = %d\n", exitflag);
disp(output.message);

Tw_d  = xsol(1:Nz);
Tw_do = xsol(Nz+1:end-1);
L1 = xsol(end);
%% ========= 3) 用解重新计算 =========
outL = sub_annular(param, Tw_d,L1);
outE = external_tube(param, Tw_do,L1);
outW = wall(param, Tw_d, Tw_do);

%% ========= 4) 结果展示 =========
z = outL.z;

figure;
plot(z, outL.T-273.15, z, outE.T_ex-273.15, z, Tw_d-273.15, z, Tw_do-273.15);
grid on
xlabel('z (m)'); ylabel('T (°C)');
legend('T_{in} (liquid)','T_{ex} (annulus)','T_{w,d}','T_{w,do}');

figure; plot(z, outL.q, z, outW.q_d); grid on
xlabel('z (m)'); ylabel('q_{inner}'); legend('q from Liquid','q_d from wall');

figure; plot(z, outE.q_ex, z, outW.q_do); grid on
xlabel('z (m)'); ylabel('q_{outer}'); legend('q_{ex} from external','q_{do} from wall');

result.param = param;
result.Tw_d  = Tw_d;
result.Tw_do = Tw_do;
result.outL  = outL;
result.outE  = outE;
result.outW  = outW;
assignin('base','HX_result',result);
