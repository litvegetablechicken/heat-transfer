function out = wall(inp)
% wall  管壁导热模型（离散形式）
% -------------------------------------------------
% Tw_d(i)  和 Tw_do(i) 与 Liquid 中的 Tw(i) 一一对应
% 逐点计算 q_d(i), q_do(i)
%
% 模型：
%   q_d *1E3*log(Di/Do)  = (k_w/(0.5*Di))*(Tw_d - Tw_do)
%   q_do*1E3*log(Di/Do) = (k_w/(0.5*Do))*(Tw_d - Tw_do)
%
% 输入 inp 结构体：
%   inp.inner_diameter      Di  [m]
%   inp.outer_diameter      Do  [m]
%   inp.thermal_cond_wall   k_w [W/(m·K)]
%   inp.Tw_d   [Nz×1]  内壁温度（离散）
%   inp.Tw_do  [Nz×1]  外壁温度（离散）
%   inp.Nz     离散点数（必须与 Tw 长度一致）
%
% 输出 out：
%   out.q_d   [Nz×1]  内侧热流密度
%   out.q_do  [Nz×1]  外侧热流密度
%   out.dT    [Nz×1]  壁温差
% -------------------------------------------------

%% -------- 1. 字段检查 --------
req = ["inner_diameter","outer_diameter","thermal_cond_wall", ...
       "Tw_d","Tw_do","Nz"];
for f = req
    if ~isfield(inp,f)
        error("wall:MissingField","inp.%s is required.",f);
    end
end

Di = mustPosScalar(inp.inner_diameter,"inner_diameter");
Do = mustPosScalar(inp.outer_diameter,"outer_diameter");
kw = mustPosScalar(inp.thermal_cond_wall,"thermal_cond_wall");
Nz = mustIntGE(inp.Nz,2,"Nz");

Tw_d  = mustVector(inp.Tw_d,"Tw_d");
Tw_do = mustVector(inp.Tw_do,"Tw_do");

if numel(Tw_d) ~= Nz || numel(Tw_do) ~= Nz
    error("wall:SizeMismatch", ...
          "Tw_d and Tw_do must have length Nz (same as Liquid).");
end

%% -------- 2. 几何检查 --------
ratio = Di/Do;
if ratio <= 0
    error("wall:InvalidGeometry","inner_diameter/outer_diameter must be > 0.");
end

logR = log(ratio);
if abs(logR) < 1e-12
    error("wall:LogZero","log(inner_diameter/outer_diameter) ≈ 0 (Di≈Do).");
end

%% -------- 3. 离散计算 --------
dT = Tw_d - Tw_do;

q_d  = (kw./(0.5*Di)) .* dT ./ (logR);
q_do = (kw./(0.5*Do)) .* dT ./ (logR);

%% -------- 4. 输出 --------
out.q_d  = q_d(:);
out.q_do = q_do(:);
out.dT   = dT(:);

out.Di = Di;
out.Do = Do;
out.kw = kw;
out.logR = logR;

end

%% ====== 工具函数 ======
function x = mustPosScalar(x,name)
if ~isscalar(x) || ~isfinite(x) || x<=0
    error("wall:InvalidInput","inp.%s must be a positive scalar.",name);
end
x = double(x);
end

function v = mustVector(v,name)
if ~isnumeric(v) || any(~isfinite(v(:)))
    error("wall:InvalidInput","inp.%s must be a finite vector.",name);
end
v = double(v(:));
end

function n = mustIntGE(n,minv,name)
if ~isscalar(n) || n~=round(n) || n<minv
    error("wall:InvalidInput","inp.%s must be integer >= %d.",name,minv);
end
n = double(n);
end
