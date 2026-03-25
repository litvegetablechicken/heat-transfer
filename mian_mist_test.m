clear
Di = 0.093;
NzM = 101;
param.geom.NzM = NzM;
param.geom.Di = Di;
param.geom.A_i  = pi*param.geom.Di^2/4;

Boundary.Mist.G_L0 = 0;
Boundary.Mist.G_V0 = 100;
Boundary.Mist.G_ED0 = 100;
param.fluid.name = 'water';

param.fluid.rhoL  = 1000;
param.fluid.rhoV  = 0.6;
param.fluid.rhoED = 1000;

param.fluid.muL = 1E-3;
param.fluid.muV = 0.015E-3;

param.fluid.kL = 0.6;
param.fluid.kV = 30E-3;

param.fluid.cpL = 4.2e3;
param.fluid.cpV = 2.0e3;
param.fluid.cpED = 4.2e3;

param.fluid.sigma = 72.9E-3;

param.fluid.DH_vap = 2257*1e3;

param.fluid.P  = 1;
param.fluid.Pc = 217.0;

param.fluid.Mw = 18.015;

param.fluid.Tbp = 100 + 273.15;
q_r = 100e3 * ones(NzM,1);   % kW/m2
x0 = 1;
opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-15, ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',2e10, ...
    'UseParallel',true);
[x_sol,~,exitflag,output] = fsolve(@(x) residual_mist(x,param, q_r, Boundary), x0, opts);
out_M = mist_q(param, q_r,Boundary,x_sol);
function F = residual_mist(x,param, q_r,Boundary)
L= x;
out_M = mist_q(param, q_r,Boundary,L);
F= out_M.G_ED(end);
end