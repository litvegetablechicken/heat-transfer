function x0 = init_double_pipe_q(param)

Nz = param.geom.Nz;

%% 温差
DT = param.external.T_in_ex - param.liquid.T_in;

%% 取液体侧物性估算 h
rho = param.liquid.density;
mu  = param.liquid.viscosity;
k   = param.liquid.thermal_cond;
cp  = param.liquid.heat_capacity;

Di = param.geom.Di;
A  = param.geom.A_i;

m  = param.liquid.m;

%% 速度与质量通量
V = m / rho;
u = V / A;
G = m / A;

%% 无量纲数
Re = rho * u * Di / mu;
Pr = mu * cp / k;

Nu = 0.023 * Re^0.8 * Pr^0.4;

%% 换热系数
h = Nu * k / Di;

%% 热流初值
q_init = h * DT;

%% 构造向量
q_d0  = q_init * ones(Nz,1);
q_do0 = q_init * ones(Nz,1);

%% 长度初值
if isfield(param.geom,'L')
    x0 = [q_d0; q_do0];
else
    L0 = [1;1];   % 默认1 m
    x0 = [q_d0; q_do0; L0];
   
end

end