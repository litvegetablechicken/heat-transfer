function out = Liquid(param, Tw_d)
Nz = param.geom.Nz;
L  = param.geom.L;

m              = param.liquid.m;
density        = param.liquid.density;
viscosity      = param.liquid.viscosity;
thermal_cond   = param.liquid.thermal_cond;
heat_capacity  = param.liquid.heat_capacity; % kJ/kgK
inner_diameter = param.geom.Di;
cross_area     = param.geom.A_i;

T_in = param.liquid.T_in;

% 网格
z  = linspace(0, L, Nz).';
dz = z(2)-z(1);

% 质量相关
V = m / density;
u = V / cross_area;
G = m / cross_area;

% 壁温向量
Tw = Tw_d(:);
if numel(Tw) ~= Nz
    error("Liquid_param:TwSize","Tw_d must be length Nz.");
end

% 无量纲数 + h（无稳定性补丁）
Re = density * u * inner_diameter / viscosity;
Pr = viscosity * (heat_capacity*1e3) / thermal_cond;
Nu = 0.023 * (Re^0.8) * (Pr^0.4);
h_const = Nu * thermal_cond / inner_diameter;

% 显式差分
T   = zeros(Nz,1);
q   = zeros(Nz,1);
h_z = zeros(Nz,1);

T(1) = T_in;
cpJ  = heat_capacity*1e3;

for k = 1:Nz-1
    h_z(k) = h_const;
    q(k)   = h_z(k) * (Tw(k) - T(k));
    dTdz   = (4*q(k)) / (G * cpJ * inner_diameter);
    T(k+1) = T(k) + dTdz*dz;
end

h_z(end) = h_const;
q(end)   = h_z(end) * (Tw(end) - T(end));

out.z  = z;
out.T  = T;
out.h  = h_z;
out.q  = q;

out.V  = V;
out.u  = u;
out.G  = G;

out.Re = Re;
out.Pr = Pr;
out.Nu = Nu;
end
