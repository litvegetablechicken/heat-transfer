clear; clc; close all

param.geom.NzSA = 101;
param.geom.Di = 0.093;
param.sub_annular.Mw = 18.015;
param.sub_annular.Pc = 217.0;
param.geom.A_i = pi*param.geom.Di^2/4;
param.sub_annular.P = 2.0;
Boundary.SA.G_L = 297;
Boundary.SA.G_V = 0;

param.sub_annular.rhoL = 1000;
param.sub_annular.rhoV = 0.6;

param.sub_annular.DH_vap = 2257e3;  % J/kg
param.sub_annular.g = 9.81;

param.liquid.Tbp = 100+273.15;
q_r = ones(param.geom.NzSA,1) * 100e3;
x0 = 1;
opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-15, ...
    'MaxIterations',20000, ...
    'MaxFunctionEvaluations',2e5, ...
    'UseParallel',false);

[x_sol, fval, exitflag, output] = fsolve(@(x) residual_SA(x, param,Boundary,q_r), x0, opts);

function F = residual_SA(x,param,Boundary,q_r)

L = x;
out = sub_annular_q(param, Boundary,q_r,L);
F = out.term;
end