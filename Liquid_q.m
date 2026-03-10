function out = Liquid_q(param, q_d, L)

Nz = param.geom.NzL;

m              = param.liquid.m;
density        = param.liquid.density;
viscosity      = param.liquid.viscosity;
thermal_cond   = param.liquid.thermal_cond;
heat_capacity  = param.liquid.heat_capacity; % J/kgK
inner_diameter = param.geom.Di;
cross_area     = param.geom.A_i;

T_in = param.liquid.T_in;

% 网格
z  = linspace(0, 1, Nz).';
dz = z(2)-z(1);

% 质量相关
V = m / density;
u = V / cross_area;
G = m / cross_area;

% 热流输入
q = q_d(:);
if numel(q) ~= Nz
    error("Liquid_param:qSize","q_d must be length Nz.");
end

% 无量纲数 + h
Re = density * u * inner_diameter / viscosity;
Pr = viscosity * (heat_capacity) / thermal_cond;
Nu = 0.023 * (Re^0.8) * (Pr^0.4);
h_const = Nu * thermal_cond / inner_diameter;

% 显式差分
T   = zeros(Nz,1);
Tw  = zeros(Nz,1);
h_z = zeros(Nz,1);

T(1) = T_in;
cpJ  = heat_capacity;

for k = 1:Nz-1

    h_z(k) = h_const;

    % 由 q 反算壁温
    Tw(k) = T(k) + q(k) / h_z(k);

    % 能量方程推进
    dTdz   = (4*L*q(k)) / (G * cpJ * inner_diameter);
    T(k+1) = T(k) + dTdz*dz;

end

h_z(end) = h_const;
Tw(end)  = T(end) + q(end) / h_z(end);

out.z  = z*L;

out.T  = T;
out.Tw = Tw;
out.h  = h_z;
out.q  = q;

out.V  = V;
out.u  = u;
out.G  = G;

out.Re = Re;
out.Pr = Pr;
out.Nu = Nu;

end