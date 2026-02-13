function out = external_tube(inp)
% external_tube  环隙/外管侧离散计算：T_ex(z), h_ex(z), q_ex(z)
% 不包含任何“稳定性补丁”：不加 eps，不加 abs(Re+1e-6)，不做 max(...)
%
% 模型：
%   m = V*density_ex
%   u*A = V
%   G = m/A
%   q_ex = h_ex*(T_ex - T_w)
%   Re = rho*u*Dh/mu
%   Pr = mu*(cp*1e3)/k
%   Nu = 0.023*(Re^0.8)*(Pr^0.4)
%   h_ex = Nu*k/Dh
%
% 能量方程：
%   dir>0:  G*cp*dTdz*(De^2-Do^2) = -(L*4*Do)*q_ex
%   dir<0: -G*cp*dTdz*(De^2-Do^2) = -(L*4*Do)*q_ex

%% 1) 字段检查（保留：防止输入缺失）
req = ["V","density_ex","viscosity_ex","thermal_cond_ex","heat_capacity_ex", ...
       "cross_area_ex","hydraulic_diameter","external_diameter","outer_diameter", ...
       "L","dir","energy_balance_cal","T_in_ex","T_w","Nz"];
for f = req
    if ~isfield(inp,f)
        error("external_tube:MissingField","inp.%s is required.",f);
    end
end

V   = mustPosScalar(inp.V,"V");
rho = mustPosScalar(inp.density_ex,"density_ex");
mu  = mustPosScalar(inp.viscosity_ex,"viscosity_ex");
k   = mustPosScalar(inp.thermal_cond_ex,"thermal_cond_ex");
cp  = mustPosScalar(inp.heat_capacity_ex,"heat_capacity_ex"); % kJ/kgK
A   = mustPosScalar(inp.cross_area_ex,"cross_area_ex");
Dh  = mustPosScalar(inp.hydraulic_diameter,"hydraulic_diameter");
De  = mustPosScalar(inp.external_diameter,"external_diameter");
Do  = mustPosScalar(inp.outer_diameter,"outer_diameter");
L   = mustPosScalar(inp.L,"L");

dir = inp.dir;
if ~isscalar(dir) || ~isfinite(dir) || dir==0
    error("external_tube:InvalidInput","inp.dir must be a finite nonzero scalar.");
end

mode = string(inp.energy_balance_cal);
Tin  = mustFiniteScalar(inp.T_in_ex,"T_in_ex");
Nz   = mustIntGE(inp.Nz,2,"Nz");

Tw_in = inp.T_w;
if isscalar(Tw_in)
    Tw = repmat(double(Tw_in), Nz, 1);
else
    Tw = double(Tw_in(:));
    if numel(Tw) ~= Nz
        error("external_tube:TwSize","inp.T_w must be scalar or length Nz.");
    end
end

if De^2 - Do^2 <= 0
    error("external_tube:InvalidGeometry","external_diameter^2 - outer_diameter^2 must be > 0.");
end

%% 2) 网格
z  = linspace(0, L, Nz).';
dz = z(2)-z(1);

%% 3) 质量相关
m_ex = V * rho;
u_ex = V / A;
G_ex = m_ex / A;

%% 4) 无量纲数 + h_ex（无稳定性补丁）
Re = rho * u_ex * Dh / mu;
Pr = mu * (cp*1e3) / k;
Nu = 0.023 * (Re^0.8) * (Pr^0.4);
h_const = Nu * k / Dh;

%% 5) 离散求解
T_ex = zeros(Nz,1);
q_ex = zeros(Nz,1);
h_ex = zeros(Nz,1);

T_ex(1) = Tin;

cpJ  = cp*1e3;
geom = (De^2 - Do^2);

% 从能量方程整理得系数：
% dir>0: dT/dz = -(L*4*Do)/(G*cp*(geom)) * q_ex
% dir<0: dT/dz = +(L*4*Do)/(G*cp*(geom)) * q_ex
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
        error("external_tube:InvalidMode", ...
              "inp.energy_balance_cal must be 'energy_bal' or 'contant_temp'/'constant_temp'.");
    end
end

% 末端点
h_ex(end) = h_const;
q_ex(end) = h_ex(end) * (T_ex(end) - Tw(end));

%% 6) 输出
out.z    = z;
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

%% ================= 输入检查工具函数（不属于数值稳定补丁） =================
function x = mustPosScalar(x,name)
if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x) || x<=0
    error("external_tube:InvalidInput","inp.%s must be a positive finite scalar.",name);
end
x = double(x);
end

function x = mustFiniteScalar(x,name)
if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x)
    error("external_tube:InvalidInput","inp.%s must be a finite scalar.",name);
end
x = double(x);
end

function n = mustIntGE(n,minv,name)
if ~isscalar(n) || ~isnumeric(n) || ~isfinite(n) || n~=round(n) || n<minv
    error("external_tube:InvalidInput","inp.%s must be integer >= %d.",name,minv);
end
n = double(n);
end
