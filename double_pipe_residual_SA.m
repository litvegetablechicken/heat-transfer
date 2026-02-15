function F = double_pipe_residual_SA(x, param)
Nz = param.geom.NzSA;

Tw_d  = x(1:Nz);
Tw_do = x(Nz+1:end-1);
L_SA = x(end);
outW = wall(param, Tw_d, Tw_do);
outL = sub_annular(param, Tw_d,L_SA);
outE = external_tube(param, Tw_do,L_SA);

r_inner = outL.q    - outW.q_d;
r_outer = outE.q_ex - outW.q_do;
r_term  =  outL.term;
F = [r_inner; r_outer;r_term];
end
