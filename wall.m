function out = wall(param, Tw_d, Tw_do)
Nz = param.geom.Nz;

Di = param.geom.Di;
Do = param.geom.Do;
kw = param.wall.kw;

Tw_d  = Tw_d(:);
Tw_do = Tw_do(:);

if numel(Tw_d) ~= Nz || numel(Tw_do) ~= Nz
    error("wall_param:SizeMismatch","Tw_d and Tw_do must be length Nz.");
end

logR = log(Di/Do);

dT   = Tw_d - Tw_do;
q_d  = (kw/(0.5*Di)) * dT / logR;
q_do = (kw/(0.5*Do)) * dT / logR;

out.q_d  = q_d;
out.q_do = q_do;
out.dT   = dT;
out.logR = logR;
end
