function out = external_tube_q(param, q_ex_in, L,Nz)
% Nz = param.geom.Nz;

V   = param.external.V;
rho = param.external.density_ex;
mu  = param.external.viscosity_ex;
k   = param.external.thermal_cond_ex;
cp  = param.external.heat_capacity_ex; % J/kgK

A   = param.geom.A_ex;
Dh  = param.geom.Dh;
De  = param.geom.De;
Do  = param.geom.Do;

dir  = param.external.dir;
mode = string(param.external.energy_balance_cal);
Tin  = param.external.T_in_ex;

% 已知热流 q_ex
q_ex = q_ex_in(:);
if numel(q_ex) ~= Nz
    error("external_tube_param:qSize","q_ex_in must be length Nz.");
end

% 网格
z  = linspace(0, 1, Nz).';
dz = z(2)-z(1);

% 质量相关
m_ex = V * rho;
u_ex = V / A;
G_ex = m_ex / A;

% 无量纲数 + h
Re = rho * u_ex * Dh / mu;
Pr = mu * (cp) / k;
Nu = 0.023 * (Re^0.8) * (Pr^0.4);
h_const = Nu * k / Dh;

% 变量初始化
T_ex = zeros(Nz,1);
Tw   = zeros(Nz,1);
h_ex = zeros(Nz,1);

T_ex(1) = Tin;

cpJ  = cp;
geom = (De^2 - Do^2);

coef = (L*4.0*Do) / (G_ex * cpJ * geom);

for i = 1:Nz-1
    h_ex(i) = h_const;

    % 已知 q_ex，反求壁温
    Tw(i) = T_ex(i) - q_ex(i) / h_ex(i);

    if mode=="energy_bal"
        if dir > 0
            dTdz = -coef * q_ex(i);
        else
            dTdz = +coef * q_ex(i);
        end

        T_ex(i+1) = T_ex(i) + dTdz*dz;

    elseif mode=="constant_temp"
        T_ex(i+1) = T_ex(i);

    else
        error("external_tube_param:InvalidMode", ...
              "energy_balance_cal must be 'energy_bal' or 'contant_temp'/'constant_temp'.");
    end
end

% 最后一个点
h_ex(end) = h_const;
Tw(end)   = T_ex(end) - q_ex(end) / h_ex(end);

out.z    = z*L;
out.T_ex = T_ex;
out.Tw   = Tw;
out.h_ex = h_ex;
out.q_ex = q_ex;

out.m_ex = m_ex;
out.u_ex = u_ex;
out.G_ex = G_ex;

out.Re = Re;
out.Pr = Pr;
out.Nu = Nu;
end