function F = double_pipe_residual(x, param)
Nz = param.geom.Nz;

Tw_d  = x(1:Nz);
Tw_do = x(Nz+1:end-1);
param.geom.L = x(end);
outW = wall(param, Tw_d, Tw_do);
outL = Liquid(param, Tw_d);
outE = external_tube(param, Tw_do);

r_inner = outL.q    - outW.q_d;
r_outer = outE.q_ex - outW.q_do;
r_Tbp = ones(param.geom.Nz,1) * (outL.T(end) - param.liquid.Tbp);
F = [r_inner; r_outer;r_Tbp];
end
