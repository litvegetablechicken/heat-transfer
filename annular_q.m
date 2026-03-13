function out = annular_q(param, q_d, L, delta_in, state0)
% annular  —— 用于替代 sub_annular 的环状流(Annular)区域模型
%
% 目标：在给定壁温 Tw_d(z) 下，沿程求解三相(液膜L/气核V/夹带滴ED)的质量通量 G，
%      并根据选定传热关联式计算壁面热流 q_r = h*(T_w - T_BP)。
%
% 说明（按你给的方程结构实现，保留方程中出现的 ABS(Re+1E-6) 等"模型自带"项）：
% - 温度：液膜 T = T_BP；气核入口 T_C(0)=T_BP
% - 质量守恒：dG_L/dz, dG_V/dz, dG_ED/dz 按给定相间源项(蒸发 V_L、夹带ED、沉积D)离散推进
% - 传热：q_r = h*(T_w - T)；h 由 h_calculation 选择（默认 correlation_0：Cooper + 简单对流）
% - 膜厚：delta_calculation 可选 'triangular_rel' 或 'vol_fraction'
% - 夹带边界：按 Barbosa-Hewitt 2002 初始夹带分数条件（可关闭）
%
% 输入：
%   param : 统一参数结构体（需包含下述字段，缺失会报错）
%   Tw_d  : Nz×1 壁面温度（内壁）
%   state0(可选): 初始状态结构体（用于指定入口 G_L0/G_V0/G_ED0, 以及开关）
%
% 输出 out：
%   out.z, out.T, out.T_C, out.q, out.h, out.delta
%   out.G_L, out.G_V, out.G_ED  (kg/m^2/s)
%   out.aux  : 其他中间量（Re_LF, Re_C, C_ED, D, ED, V_L 等）
%
% -------------------------------------------------------------

%% ===================== 0) 基本取参 =====================
Nz = param.geom.Nz;

Di = param.geom.Di;
Area  = param.geom.A_i;

Tw = Tw_d(:);
if numel(Tw) ~= Nz
    error("annular:TwSize","Tw_d must be length Nz.");
end


z = linspace(0,1,Nz).';
dz    = z(2)-z(1);

% 相：L, V, ED
rhoL = param.annular.rhoL;   muL = param.annular.muL;   kL  = param.annular.kL;   cpL = param.annular.cpL; % cp: kJ/kgK
rhoV = param.annular.rhoV;   muV = param.annular.muV;   kV  = param.annular.kV;   cpV = param.annular.cpV;
rhoED= param.annular.rhoED;  cpED= param.annular.cpL;  % ED 常取与液相相同也可，但这里显式给

sigma = param.annular.surface_tension;     % N/m
DHvap = param.annular.DH_vap;              % kJ/kg
g     = param.annular.g;

% 饱和/Cooper 需要的参数
P   = param.annular.P;      % bar 或者与 Pc 同单位（只要 P/Pc 无量纲）
Pc  = param.annular.Pc;
Mw  = param.annular.Mw;     % kg/kmol 或 g/mol 只要一致（指数用）
% Antoine 常数如需 dPsatdT 可再加，这里只实现 correlation_0 主流程，不强制 Antoine

% 设置默认计算选择开关
if isfield(param,'annular') && isfield(param.annular,'delta_calculation')
    delta_calculation = string(param.annular.delta_calculation);
else
    delta_calculation = "triangular_rel";
end

if isfield(param,'annular') && isfield(param.annular,'entraintment_calculation')
    entraintment_calculation = string(param.annular.entraintment_calculation);
else
    entraintment_calculation = "cal_entraintment"; 
end

if isfield(param,'annular') && isfield(param.annular,'h_calculation')
    h_calculation = string(param.annular.h_calculation);
else
    h_calculation = "correlation_0";
end

% 入口边界/初值
if nargin < 5 || isempty(state0)
    state0 = struct();
end

if ~isfield(state0,'calculate_entraintment_0'); state0.calculate_entraintment_0 = true; end

G_L  = zeros(Nz,1);  G_V  = zeros(Nz,1);  G_ED = zeros(Nz,1);

G_L(1) = Boundary.A.G_L;
G_V(1) = Boundary.A.G_V;
% ---- BOUNDARY：G_ED(0) ----
if G_V(1)==0 || ~state0.calculate_entraintment_0
    G_ED(1) = 1e-10;
else
    % Barbosa-Hewitt 2002 初始夹带分数条件：
    % (1/(0.95E-2 + 342.55E-2*sqrt((rhoL*G_L)/(rhoV*G_V))*Di^2) - 1)*G_ED = G_L
    coefBH = (1.0/(0.95e-2 + 342.55e-2 * sqrt((rhoL*G_L(1))/(rhoV*G_V(1))) * (Di^2)) - 1.0);
    % coefBH * G_ED = G_L  ->  G_ED = G_L/coefBH
    G_ED(1) = G_L(1) / coefBH;
end

% 温度：T = T_BP；T_C(0)=T_BP
T_BP = param.liquid.Tbp;     % 你主程序里 param.liquid.Tbp
T    = ones(Nz,1)*T_BP;
T_C  = zeros(Nz,1); T_C(1) = T_BP;

% 输出：热流与传热系数
q_r  = zeros(Nz,1);
h    = zeros(Nz,1);
if isscalar(delta_in)
    delta = ones(Nz,1) * delta_in;
else
    delta = delta_in(:);
    if numel(delta) ~= Nz
        error("annular:DeltaSize","delta_in must be scalar or length Nz.");
    end
end

% 一些中间量
D   = zeros(Nz,1);
ED  = zeros(Nz,1);
V_L = zeros(Nz,1);
C_ED= zeros(Nz,1);
Re_C= zeros(Nz,1);
Re_LF = zeros(Nz,1);
Pr_LF = zeros(Nz,1);
delta_new = zeros(Nz,1);


%% ===================== 1) 主循环：沿程推进 =====================
for k = 1:Nz

    % ---------- (A) 相速度/体积流量/质量分数（按你定义） ----------
    % m = G*A  (每相)
    mL  = G_L(k)*Area;
    mV  = G_V(k)*Area;
    mED = G_ED(k)*Area;

    % 体积流量 V = m/rho
    VL  = mL / rhoL;
    VV  = mV / rhoV;
    VED = mED/ rhoED;

    % 速度 u = V/A
    uL  = VL / Area;
    uV  = VV / Area;
    uED = VED/ Area;

    % 体积分数 w_frac： w_frac(ph)*SUM(m/rho) = m(ph)/rho(ph)
    denomVol = (mL/rhoL + mV/rhoV + mED/rhoED);
    wL = (mL/rhoL) / denomVol;
    wV = (mV/rhoV) / denomVol;
    wED= (mED/rhoED)/ denomVol;

    % 质量分数 v_frac： v_frac(ph)*SUM(m) = m(ph)
    denomM = (mL + mV + mED);
    vL = mL/denomM; vV = mV/denomM; vED = mED/denomM;

    % ---------- (B) 夹带分数 E_frac*(mED+mL)=mED ----------
    E_frac = mED / (mED + mL);

    % ---------- (C) 液膜厚度 delta ----------
    if delta_calculation == "vol_fraction"
        % delta = 0.5*Di*w_frac(L)
        delta(k) = 0.5*Di*wL;
    else
        % triangular_rel：需要 tau_I, f_I, f_LF, Re_C, Re_LF 等
        % 气核混合密度 rho_C：rho_C*(mED/rhoED + mV/rhoV) = mED + mV
        rhoC = (mED + mV) / (mED/rhoED + mV/rhoV);

        % u_C*((1 - 2*delta/Di)^2) = G_ED/rhoED + G_V/rhoV
        % -> u_C = (G_ED/rhoED + G_V/rhoV) / (1-2*delta/Di)^2
        % 这里 delta 与 u_C 耦合，用小迭代求 delta
        delta_it = delta(max(k-1,1));
        uC = (G_ED(k)/rhoED + G_V(k)/rhoV) / ((1 - 2*delta_it/Di)^2);

        % Re_C*muV = rhoC*uC*Di
        ReC_it = (rhoC*uC*Di)/muV;

        % f_SC = 0.079/((ABS(Re_C+1E-6)^0.25))
        fSC = 0.079/(abs(ReC_it + 1e-6)^0.25);

        % f_I = f_SC*(1 + 24*(delta/Di)*(rhoL/rhoV)^(1/3))
        fI = fSC*(1 + 24*(delta_it/Di)*(rhoL/rhoV)^(1/3));

        % tau_I = 0.5*rhoC*f_I*(u_C - u_L)^2
        tauI = 0.5*rhoC*fI*(uC - uL)^2;

        % Re_LF*muL = rhoL*uL*(4*delta)
        ReLF_it = (rhoL*uL*(4*delta_it))/muL;

        % f_LF = 0.079/((ABS(Re_LF+1E-6)^0.25))
        fLF = 0.079/(abs(ReLF_it + 1e-6)^0.25);

        % (G_L^2)*f_LF = ((4*delta/Di)^2)*(2*tau_I*rhoL)
        % -> delta_new = (Di/4)*sqrt( (G_L^2*f_LF)/(2*tau_I*rhoL) )
        delta_new(k) = (Di/4)*sqrt( (G_L(k)^2 * fLF) / (2*tauI*rhoL) );


        % 把一些中间量存一下
        Re_C(k)  = (rhoC*((G_ED(k)/rhoED + G_V(k)/rhoV)/((1 - 2*delta_new(k)/Di)^2))*Di)/muV;
        Re_LF(k) = (rhoL*uL*(4*delta_new(k)))/muL;
        Pr_LF(k) = muL*(cpL)/kL;
    end

    % ---------- (D) 传热系数 h 与热流 q_r：按你的选择器 ----------
    DT = (Tw(k) - T(k));
    h_nb = (55.0*((P/Pc)^0.12)*(-log10(P/Pc))^(-0.55)*(Mw^(-0.5))*((q_d)^0.67)); % W/m2K

    % 对流项：Nu_LF = 0.0133*(Re_LF^0.69)*(Pr_LF^0.4)
    % h_conv*1E3 与 Nu_LF*thermal_cond("L")=(h_conv*1E3)*delta  => h_conv = Nu_LF*kL/(1e3*delta)
    if delta_new(k) < 0
        error("annular:deltaNonPositive","delta <= 0 at k=%d",k);
    end
    ReLF = (rhoL*uL*(4*delta_new(k)))/muL;
    PrLF = muL*(cpL)/kL;
    NuLF = 0.0133*(ReLF^0.69)*(PrLF^0.4);
    if delta_new(k) ~= 0
        h_conv = (NuLF*kL)/(delta_new(k)); % W/m2K（注意这里严格按你式子里的 1E3）
    else
        h_conv = 0;
    end
    % 选择器（这里实现你最常用的 correlation_0 / simple；其他可按需扩展）
    switch h_calculation
        case "simple"
            h_loc = 300.0; % W/m2K（仅初始化用途）
        case "correlation_0"
            % h^3 = h_nb^3 + h_conv^3
            h_loc = (h_nb^3 + h_conv^3)^(1/3);
        otherwise
            error("annular:UnsupportedH","h_calculation=%s not implemented in this MATLAB stub.",h_calculation);
    end

    q_r(k) = h_loc * DT;
    h(k) = h_loc;

    % ---------- (E) 蒸发通量 V_L = q_r/DH_vap ----------
    % 你写的是：V_L = q_r/DH_vap；q_CS = V_L*DH_vap；因此 DH_vap 用 kJ/kg 时要注意单位
    % q_r [W/m2] = J/s/m2, DH_vap [kJ/kg] => 1 kJ/kg = 1000 J/kg
    V_L(k) = q_r(k) / (DHvap);  % kg/(m2·s)

    % ---------- (F) 夹带/沉积（总开关） ----------
    % Droplet concentration： C_ED*INTEGRAL(i in [V,ED]; m(i)/rho(i)) = mED
    denomCoreVol = (mV/rhoV + mED/rhoED);
    C_ED(k) = mED / denomCoreVol; % kg/m3（按原式）

    % Deposition k_D 条件
    if (C_ED(k)/rhoV) <= 0.3
        kD = 0.18*sqrt(sigma/(rhoV*Di));
    else
        kD = 0.083*((rhoV/C_ED(k))^0.65)*sqrt(sigma/(rhoV*Di));
    end

    % Entraintment flux ED_cal 条件
    GC = (muL/Di)*exp(5.8504 + 0.4249*(muV/muL)*sqrt(rhoL/rhoV));
    if G_L(k) > GC
        ED_cal = 5.75e-5*G_V(k) * ( ( ((rhoL*Di)/(sigma*(rhoV^2))) * ((G_L(k)-GC)^2) )^0.316 );
    else
        ED_cal = 0.0;
    end

    % Selector
    switch entraintment_calculation
        case "no_entraintment"
            D(k)  = 0.0;
            ED(k) = 0.0;
        case "cal_entraintment"
            D(k)  = kD*C_ED(k);  % deposition flux
            ED(k) = ED_cal;      % entrainment flux
        otherwise
            error("annular:UnsupportedEntrain","entraintment_calculation=%s not supported.",entraintment_calculation);
    end

    % ---------- (G) 质量守恒：沿 z 推进（显式一阶） ----------
    % 你给的连续方程： PARTIAL(G, z)*Di = L*4*(source)
    % 注意：你原方程使用 i∈[0,1] 的归一化坐标，左侧是 ∂G/∂z * Di
    % 我们这里 z 已经是 Lz，dz 是实际长度增量，所以直接用：
    % dG/dz = (L*4/Di)*source
    if k < Nz
        % V_ED=0（按你后面写的 V_ED = 0.0）
        V_ED = 0.0;

        % Liquid: dG_L/dz*Di = L*4*(D - ED - V_L)
        % 你还有 m=0 的 MAX(0,...) 条件，这里按式实现
        srcL = (D(k) - ED(k) - V_L(k));
        if mL == 0.0
            srcL = max(0.0, srcL);
        end
        dG_L_dz = (L*4.0/Di) * srcL;

        % Vapor: dG_V/dz*Di = L*4*(V_L + V_ED)
        dG_V_dz = (L*4.0/Di) * (V_L(k) + V_ED);

        % Droplets: dG_ED/dz*Di = L*4*(ED - D - V_ED)
        srcED = (ED(k) - D(k) - V_ED);
        if mED == 0.0
            srcED = max(0.0, srcED);
        end
        dG_ED_dz = (L*4.0/Di) * srcED;

        G_L(k+1)  = G_L(k)  + dG_L_dz  * dz;
        G_V(k+1)  = G_V(k)  + dG_V_dz  * dz;
        G_ED(k+1) = G_ED(k) + dG_ED_dz * dz;
    
        % 气核能量方程这里先按 equilibrium（dT_C/dz=0）实现：T_C=T_BP
        T_C(k+1) = T_BP;
    end
end

%% ===================== 2) 输出 =====================
out.z    = z*L;
out.T    = T;          % 液膜温度（T_BP）
out.T_C  = T_C;        % 气核温度（此实现为 T_BP 常数）
out.q    = q_r;        % 用于与 wall.q_d 做残差匹配（替代 Liquid.out.q）
out.h    = h;
out.delta= delta;
out.delta_term = (delta - delta_new)./delta;
out.term = G_L(end);
out.G_L  = G_L;
out.G_V  = G_V;
out.G_ED = G_ED;

out.aux.v_frac = [ ];  %  % 如需可扩展输出 v_frac/w_frac
out.aux.w_frac = [ ];
out.aux.D      = D;
out.aux.ED     = ED;
out.aux.V_L    = V_L;
out.aux.C_ED   = C_ED;
out.aux.Re_C   = Re_C;
out.aux.Re_LF  = Re_LF;
out.aux.Pr_LF  = Pr_LF;

end
