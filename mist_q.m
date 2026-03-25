function out = mist_q(param, q_r,Boundary,L)

Nz = param.geom.NzM;
z  = linspace(0,1,Nz);
dz = z(2)-z(1);

Di = param.geom.Di;

% ---------------------------
% 物性
% ---------------------------
rhoL = param.fluid.rhoL;
rhoV = param.fluid.rhoV;
muL  = param.fluid.muL;
muV  = param.fluid.muV;
kV   = param.fluid.kV;
cpV  = param.fluid.cpV;
DH   = param.fluid.DH_vap;
T = param.fluid.Tbp;
A = param.geom.A_i;   % 截面积

% ---------------------------
% 初始化变量
% ---------------------------
G_L  = zeros(Nz,1);
G_V  = zeros(Nz,1);
G_ED = zeros(Nz,1);

% 入口条件
G_L(1)  = Boundary.Mist.G_L0;
G_V(1)  =  Boundary.Mist.G_V0;
G_ED(1) =  Boundary.Mist.G_ED0;

vL = zeros(Nz,1);
vV = zeros(Nz,1);
vED = zeros(Nz,1);
V_ED = zeros(Nz,1);
Tw = zeros(Nz,1);
% ---------------------------
% 蒸发项（常数）
% ---------------------------

% ---------------------------
% 主循环（有限差分）
% ---------------------------
for i = 1:Nz-1
    V_ED(i) = q_r(i) / DH;

    % ====== Step 1: 更新质量通量 ======
    dG_L  = 0;
    dG_V  = 4*L/Di * (V_ED(i));
    dG_ED = 4*L/Di * (-V_ED(i));

    G_L(i+1)  = G_L(i)  + dz * dG_L;
    G_V(i+1)  = G_V(i)  + dz * dG_V;
    G_ED(i+1) = G_ED(i) + dz * dG_ED;

    % ====== Step 2: 质量 → 体积分数 ======
    mL  = G_L(i)*A;
    mV  = G_V(i)*A;
    mED = G_ED(i)*A;

    m_tot = mL + mV + mED;

    % v_frac
    vL(i)  = mL  / m_tot;
    vV(i)  = mV  / m_tot;
    vED(i) = mED / m_tot;
    if i ==Nz-1
        mL  = G_L(i+1)*A;
        mV  = G_V(i+1)*A;
        mED = G_ED(i)*A;
        vL(i+1)  = mL  / m_tot;
        vV(i+1)  = mV  / m_tot;
        vED(i+1) = mED / m_tot;
    end

end


% ---------------------------
% 传热计算（Nu, Re, Pr）
% ---------------------------
Re = zeros(Nz,1);
Pr = zeros(Nz,1);
Nu = zeros(Nz,1);
h  = zeros(Nz,1);

for i = 1:Nz

    Re(i) = ( (G_V(i)+G_ED(i)) * Di * (rhoL * vV(i) + rhoV * (1 - vV(i)))) / (muV * rhoL);

    Pr(i) = muV * cpV / kV;

    Nu(i) = 0.023 * Re(i)^0.8 * Pr(i)^0.4;

    h(i)  = Nu(i)*kV/Di;
    Tw(i) = q_r(i)  / h(i) + T;
end

% ---------------------------
% 输出
% ---------------------------
out.z   = z * L;
out.G_L = G_L;
out.G_V = G_V;
out.G_ED = G_ED;
out.h = h;
out.Re = Re;
out.Nu = Nu;
out.Tw = Tw;
end