function out = sub_annular_q(param,Boundary,q_r,L)
% sub_annaular  (入口纯液版，无 x_in)

%% ---- 0) 读参数 ----
Nz = param.geom.NzSA;
Di = param.geom.Di;
A  = param.geom.A_i;

if numel(q_r) ~= Nz
    error("sub_annaular:TwSize","Tw_d must be length Nz.");
end


Tbp = param.liquid.Tbp;
T   = ones(Nz,1) * Tbp;


P    = param.sub_annular.P;
Pc   = param.sub_annular.Pc;
Mw   = param.sub_annular.Mw;
DH   = param.sub_annular.DH_vap;  % J/kg
g    = param.sub_annular.g;
rhoL = param.sub_annular.rhoL;
rhoV = param.sub_annular.rhoV;

%% ---- 1) 网格 ----
z  = linspace(0, 1, Nz).';
dz = z(2) - z(1);

%% ---- 2) 入口条件：纯液 ----
G_L = zeros(Nz,1);
G_V = zeros(Nz,1);
G_L(1) = Boundary.SA.G_L;
G_V(1) = Boundary.SA.G_V;

%% ---- 3) 热流与换热系数（Cooper 核态沸腾）----
% q_r = h*(Tw - Tbp)
% h*1E3 = 55*(P/Pc)^0.12*(-log10(P/Pc))^(-0.55)*Mw^(-0.5)*((q_r*1E3)^0.67)
Tw = zeros(Nz,1);   % K
C = 55.0 * ( (P/Pc)^0.12 ) * ( -log10(P/Pc) )^(-0.55) * ( Mw^(-0.5) );

h   = zeros(Nz,1);   % W/(m^2 K)
V_L = zeros(Nz,1);   % kg/(m^2 s) 这里按"蒸发质量通量"理解（DH为J/kg）

for k = 1:Nz

    h(k)   = C * q_r(k)^0.67;
    V_L(k) = q_r(k) / DH;
    Tw(k) = q_r(k) / h(k) + T(k);

end

%% ---- 4) 两相质量守恒（离散积分）----
% dG_L/dz * Di = L*4*(-V_L)
% dG_V/dz * Di = L*4*(+V_L)
for k = 1:Nz-1
    dG_L_dz = (L*4.0*(-V_L(k))) / Di;
    dG_V_dz = (L*4.0*(+V_L(k))) / Di;

    G_L(k+1) = G_L(k) + dG_L_dz * dz;
    G_V(k+1) = G_V(k) + dG_V_dz * dz;
end

%% ---- 5) 速度与无量纲气相速度 ----
u_L = G_L / rhoL;
u_V = G_V / rhoV;

u_ast_C = u_V * sqrt( rhoV / ( g*Di*(rhoL - rhoV) ) );

%% ---- 6) 体积分数/质量分数（按你给的关系）----
G_sum = G_L + G_V;
V_sum = (G_L/rhoL) + (G_V/rhoV);

v_frac_L = G_L ./ G_sum;
v_frac_V = G_V ./ G_sum;

w_frac_L = (G_L/rhoL) ./ V_sum;
w_frac_V = (G_V/rhoV) ./ V_sum;

%% ---- 7) 输出----
out.z = z*L;
out.T = T;

out.h = h;
out.Tw = Tw;   
out.V_L = V_L;

out.G_L = G_L;
out.G_V = G_V;

out.u_L = u_L;
out.u_V = u_V;
out.u_ast_C = u_ast_C;

out.v_frac_L = v_frac_L;
out.v_frac_V = v_frac_V;
out.w_frac_L = w_frac_L;
out.w_frac_V = w_frac_V;

out.term = u_ast_C(end) - 1;  % 终止条件残差

end
