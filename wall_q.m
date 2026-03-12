function out = wall_q(param, q_d, q_do)

Di = param.geom.Di;
Do = param.geom.Do;
kw = param.wall.kw;

q_d  = q_d(:);
q_do = q_do(:);


logR = log(Di/Do);

% 由内侧热流计算温差
dT_from_qd  = q_d  * logR * (0.5*Di) / kw;

% 由外侧热流计算温差
dT_from_qdo = q_do * logR * (0.5*Do) / kw;

% 平均温差（作为最终壁温差）
dT = 0.5 * (dT_from_qd + dT_from_qdo);

% 两种计算方式的差值（用于检查一致性）
dT_diff = dT_from_qd - dT_from_qdo;

out.q_d  = q_d;
out.q_do = q_do;
out.dT = dT;
out.dT_diff = dT_diff;

end