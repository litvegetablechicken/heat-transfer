function out = external_tube(param, Tw_do,L)
Nz = param.geom.Nz;


V   = param.external.V;
rho = param.external.density_ex;
mu  = param.external.viscosity_ex;
k   = param.external.thermal_cond_ex;
cp  = param.external.heat_capacity_ex; % kJ/kgK

A   = param.geom.A_ex;
Dh  = param.geom.Dh;
De  = param.geom.De;
Do  = param.geom.Do;

dir  = param.external.dir;
mode = string(param.external.energy_balance_cal);
Tin  = param.external.T_in_ex;

% 壁温
Tw = Tw_do(:);
if numel(Tw) ~= Nz
    error("external_tube_param:TwSize","Tw_do must be length Nz.");
end

% 网格
z  = linspace(0, 1, Nz).';
dz = z(2)-z(1);

% 质量相关
m_ex = V * rho;
u_ex = V / A;
G_ex = m_ex / A;

% 无量纲数 + h（无稳定性补丁）
Re = rho * u_ex * Dh / mu;
Pr = mu * (cp*1e3) / k;
Nu = 0.023 * (Re^0.8) * (Pr^0.4);
h_const = Nu * k / Dh;

% 离散
T_ex = zeros(Nz,1);
q_ex = zeros(Nz,1);
h_ex = zeros(Nz,1);

T_ex(1) = Tin;

cpJ  = cp*1e3;
geom = (De^2 - Do^2);

coef = (L*4.0*Do) / (G_ex * cpJ * geom);

for i = 1:Nz-1
    h_ex(i) = h_const;

    if mode=="energy_bal"
        q_ex(i) = h_ex(i) * (T_ex(i) - Tw(i));

        if dir > 0
            dTdz = -coef * q_ex(i);
        else
            dTdz = +coef * q_ex(i);
        end

        T_ex(i+1) = T_ex(i) + dTdz*dz;

    elseif mode=="contant_temp" || mode=="constant_temp"
        q_ex(i)   = h_ex(i) * (T_ex(i) - Tw(i));
        T_ex(i+1) = T_ex(i);

    else
        error("external_tube_param:InvalidMode", ...
              "energy_balance_cal must be 'energy_bal' or 'contant_temp'/'constant_temp'.");
    end
end

h_ex(end) = h_const;
q_ex(end) = h_ex(end) * (T_ex(end) - Tw(end));

out.z    = z*L;
out.T_ex = T_ex;
out.h_ex = h_ex;
out.q_ex = q_ex;

out.m_ex = m_ex;
out.u_ex = u_ex;
out.G_ex = G_ex;

out.Re = Re;
out.Pr = Pr;
out.Nu = Nu;
end
