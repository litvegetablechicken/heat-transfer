function x0 = init_double_pipe_q(param,Th,Tl,zone)

if zone =='L'
    Nz = param.geom.NzL;
end
if zone == 'SA'
    Nz = param.geom.NzSA;
end
if zone == 'A'
    Nz = param.geom.NzA;

end
if zone =='M'
    Nz = param.geom.NzM;
end
if zone =='V'
    Nz = param.geom.NzV;
end
%% 温差
DT = Th-Tl;

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
delta0 = [];
%% 长度初值

if zone =='L'
    L = 1;   % 默认1 m
end
if zone == 'SA'
    L = 1;   % 默认1 m
end
if zone == 'A'
    % q_d0 = param.annular.q_d0;
    % q_do0 = param.annular.q_do0;
    % delta0 = param.annular.delta0;
    A = param.fluid.DH_vap * m / h /20;
    L = A / pi /Di;   % 默认1 m
end
if zone =='M'
    L = 1;   % 默认1 m
end
if zone =='V'
    L = [];   % 默认1 m
end
x0 = [q_d0; q_do0;delta0;L];

end