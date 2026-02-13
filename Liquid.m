function out = Liquid(inp)
% Liquid  内管液相沿程离散计算：T(z), h(z), q(z)
% 输入：inp 结构体
% 输出：out.z, out.T, out.h, out.q 以及 Re/Pr/Nu/u/G 等

%% -------------------------
% 1) 必需字段检查
%% -------------------------
req = ["m","density","viscosity","thermal_cond","heat_capacity", ...
       "cross_area","inner_diameter","L","T_in","T_w","Nz"];
for f = req
    if ~isfield(inp, f)
        error("Liquid:MissingField", "inp.%s is required.", f);
    end
end

% 基本合法性检查（标量/正数/有限/整数等）
m              = mustScalarPosFinite(inp.m,             "m");
density        = mustScalarPosFinite(inp.density,       "density");
viscosity      = mustScalarPosFinite(inp.viscosity,     "viscosity");
thermal_cond   = mustScalarPosFinite(inp.thermal_cond,  "thermal_cond");
heat_capacity  = mustScalarPosFinite(inp.heat_capacity, "heat_capacity"); % kJ/kgK
cross_area     = mustScalarPosFinite(inp.cross_area,    "cross_area");
inner_diameter = mustScalarPosFinite(inp.inner_diameter,"inner_diameter");
L              = mustScalarPosFinite(inp.L,             "L");

T_in = mustScalarFinite(inp.T_in, "T_in");
T_w  = mustFinite(inp.T_w, "T_w");

Nz = mustScalarIntGE(inp.Nz, 2, "Nz");

%% -------------------------
% 2) 网格
%% -------------------------
z  = linspace(0, L, Nz).';
dz = z(2) - z(1);

%% -------------------------
% 3) 质量相关量
%% -------------------------
V = m / density;         % m^3/s
u = V / cross_area;      % m/s
G = m / cross_area;      % kg/(m^2·s)

%% -------------------------
% 4) 壁温处理成向量
%% -------------------------
if isscalar(T_w)
    Tw = repmat(T_w, Nz, 1);
else
    Tw = T_w(:);
    if numel(Tw) ~= Nz
        error("Liquid:TwSize", "inp.T_w must be scalar or a vector with length Nz.");
    end
end

%% -------------------------
% 5) 无量纲数 + 换热系数（Dittus-Boelter）
%    （已去除 max(...,eps) 的“稳定性补丁”）
%% -------------------------
Re = density * u * inner_diameter / viscosity;
Pr = viscosity * (heat_capacity*1e3) / thermal_cond;
Nu = 0.023 * (Re^0.8) * (Pr^0.4);

h_const = Nu * thermal_cond / inner_diameter; % W/(m^2·K)

%% -------------------------
% 6) 离散沿程求解（显式差分）
%% -------------------------
T   = zeros(Nz,1);
q   = zeros(Nz,1);
h_z = zeros(Nz,1);

T(1) = T_in;
cpJ  = heat_capacity * 1e3;

for k = 1:Nz-1
    h_z(k) = h_const;
    q(k)   = h_z(k) * (Tw(k) - T(k));  % W/m^2

    % 严格按你给的式子：G*cp*dT/dz*Di = 4*L*q
    dTdz = (4*q(k)) / (G * cpJ * inner_diameter);
    T(k+1) = T(k) + dTdz * dz;
end

h_z(end) = h_const;
q(end)   = h_z(end) * (Tw(end) - T(end));

%% -------------------------
% 7) 输出
%% -------------------------
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

%% =========================
% 辅助校验函数（输入检查，不属于稳定性补丁）
%% =========================
function x = mustScalarPosFinite(x, name)
    if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x) || x <= 0
        error("Liquid:InvalidInput", "inp.%s must be a positive finite scalar.", name);
    end
    x = double(x);
end

function x = mustScalarFinite(x, name)
    if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x)
        error("Liquid:InvalidInput", "inp.%s must be a finite scalar.", name);
    end
    x = double(x);
end

function x = mustFinite(x, name)
    if ~isnumeric(x) || any(~isfinite(x(:)))
        error("Liquid:InvalidInput", "inp.%s must contain only finite numeric values.", name);
    end
    x = double(x);
end

function n = mustScalarIntGE(n, minVal, name)
    if ~isscalar(n) || ~isnumeric(n) || ~isfinite(n) || n ~= round(n) || n < minVal
        error("Liquid:InvalidInput", "inp.%s must be an integer scalar >= %d.", name, minVal);
    end
    n = double(n);
end
