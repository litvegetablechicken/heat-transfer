clear
clc

%% =========================================================
% 1  Geometry
%% =========================================================

param.geom.NzA = 101;      % grid number
param.geom.Di  = 0.0093;   % inner diameter [m]

% L = 1.0;                   % tube length [m]

%% =========================================================
% 2  Fluid properties (example: R245fa)
%% =========================================================

param.fluid.name = 'water';

param.fluid.rhoL  = 1000;
param.fluid.rhoV  = 1000;
param.fluid.rhoED = 0.6;

param.fluid.muL = 1E-3;
param.fluid.muV = 0.015E-3;

param.fluid.kL = 0.6;
param.fluid.kV = 30E-3;

param.fluid.cpL = 4.2e3;
param.fluid.cpV = 2.0e3;
param.fluid.cpED = 4.2e3;

param.fluid.sigma = 72.9E-3;

param.fluid.DHvap = 2257;

param.fluid.P  = 101.325;
param.fluid.Pc = 217.0e5;

param.fluid.Mw = 18.015;

param.fluid.Tbp = 100 + 273.15;

% Antoine coefficients
param.fluid.A_ant = 18.015;
param.fluid.B_ant = 1810.94;
param.fluid.C_ant = -28.665;

%% =========================================================
% 3  Model options
%% =========================================================

param.model.entraintment_calculation = "cal_entraintment";

param.model.h_calculation = "correlation_0";

%% =========================================================
% 4  Initial state
%% =========================================================

param.state0.G_L0 = 297;

param.state0.G_V0 = 50;

param.state0.calculate_entraintment_0 = true;

%% =========================================================
% 5  Boundary conditions
%% =========================================================

Boundary.A.G_L = param.state0.G_L0;

Boundary.A.G_V = param.state0.G_V0;

%% =========================================================
% 6  Wall heat flux initial guess
%% =========================================================

Nz = param.geom.NzA;

q_r = 1000 * ones(Nz,1);   % kW/m2

%% =========================================================
% 7  Run annular model
%% =========================================================
x0 = 1;
opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-15, ...
    'MaxIterations',20000, ...
    'MaxFunctionEvaluations',2e5, ...
    'UseParallel',false);
[x_sol, fval, exitflag, output] = fsolve(@(x) A_residual(x,param,q_r,Boundary), x0, opts);
%% =========================================================
% 8  Plot results
%% =========================================================
L = x_sol;
out = annular_heat_transfer(param, q_r, Boundary, L);
figure
plot(out.z,out.G_L,'LineWidth',2)
hold on
plot(out.z,out.G_V,'LineWidth',2)
plot(out.z,out.G_ED,'LineWidth',2)
legend('G_L','G_V','G_{ED}')
xlabel('z (m)')
ylabel('Mass flux (kg/m^2s)')
grid on


figure
plot(out.z,out.delta,'LineWidth',2)
xlabel('z (m)')
ylabel('\delta (m)')
grid on


figure
plot(out.z,out.T_w,'LineWidth',2)
xlabel('z (m)')
ylabel('Wall temperature (K)')
grid on


figure
plot(out.z,out.h,'LineWidth',2)
xlabel('z (m)')
ylabel('h (kW/m^2K)')
grid on

disp("Total heat transfer:")
disp(out.Q)
function F = A_residual(x,param,q_r,Boundary)
L = x;
out = annular_heat_transfer(param, q_r, Boundary, L);
F = out.term;
end