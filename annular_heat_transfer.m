function out = annular_heat_transfer(param, T_w, delta)
%ANNULAR_HEAT_TRANSFER
% Grid-based MATLAB implementation of annular-region equations.
%
% INPUT
%   param : struct, must contain geometry / fluid / model / initial fields
%   T_w   : wall temperature vector [K], size Nz x 1
%   delta : liquid film thickness vector [m], size Nz x 1  <-- modified
%
% OUTPUT
%   out   : struct
%       .q        local heat flux [kW/m^2]
%       .Q        total heat transfer [W]
%       .h        local HTC [kW/m^2-K]
%       .delta    liquid film thickness [m]
%       .G_L      liquid mass flux [kg/m^2/s]
%       .G_V      vapor mass flux [kg/m^2/s]
%       .G_ED     entrained droplet mass flux [kg/m^2/s]
%       .T        liquid film temperature [K]
%       .T_C      core temperature [K]
%       .z        axial coordinate [m]

%% -----------------------------
% 0) Basic checks and extraction
% ------------------------------
Nz = param.geom.Nz;
L  = param.L;
Di = param.geom.Di;

if isrow(T_w),   T_w   = T_w(:);   end
if isrow(delta), delta = delta(:); end

assert(length(T_w) == Nz,   'T_w must have length Nz.');
assert(length(delta) == Nz, 'delta must have length Nz.');

z  = linspace(0, 1, Nz).';         % normalized coordinate
dz = z(2) - z(1);
A  = pi*(Di/2)^2;                  %#ok<NASGU>

% Shorthand
f = param.fluid;      % fluid properties / constants
m = param.model;      % model switches
s = param.state0;     % initial state

% Required fluid properties
rhoL = f.rhoL; rhoV = f.rhoV; rhoED = f.rhoED;
muL  = f.muL;  muV  = f.muV;
kL   = f.kL;   kV   = f.kV;
cpL  = f.cpL;  cpV  = f.cpV; cpED = f.cpED;    % [kJ/kg-K] recommended
sigma = f.sigma;
DHvap = f.DHvap;      % [kJ/kg]
P     = f.P;          % consistent with your correlation
Pc    = f.Pc;
Mw    = f.Mw;

% Antoine constants if needed by correlation_2 / dPsatdT
A_ant = f.A_ant;
B_ant = f.B_ant;
C_ant = f.C_ant;

% boiling point / saturation temperature profile
if isscalar(f.Tbp)
    T_BP = f.Tbp * ones(Nz,1);
else
    T_BP = f.Tbp(:);
end

g = 9.81;

%% -----------------------------
% 1) Allocate variables
% ------------------------------
G_L  = zeros(Nz,1);   % kg/m^2/s
G_V  = zeros(Nz,1);
G_ED = zeros(Nz,1);

u_L  = zeros(Nz,1);
u_V  = zeros(Nz,1);
u_ED = zeros(Nz,1);
m_L = zeros(Nz,1);
m_V = zeros(Nz,1);
m_ED = zeros(Nz,1);
T     = T_BP;         % equilibrium model: liquid film temperature = T_BP
T_C   = zeros(Nz,1);  % gas core temperature

rho_C = zeros(Nz,1);
u_C   = zeros(Nz,1);
Re_C  = zeros(Nz,1);
Re_G  = zeros(Nz,1);
We_G  = zeros(Nz,1);
tau_I = zeros(Nz,1);
f_SC  = zeros(Nz,1);
f_I   = zeros(Nz,1);
f_LF  = zeros(Nz,1);
Re_LF = zeros(Nz,1);

C_ED  = zeros(Nz,1);
k_D   = zeros(Nz,1);
G_C   = zeros(Nz,1);
ED_cal = zeros(Nz,1);

D   = zeros(Nz,1);
ED  = zeros(Nz,1);
V_L = zeros(Nz,1);
V_ED = zeros(Nz,1);

X_vap = zeros(Nz,1);
X_tt  = zeros(Nz,1);

Pr_C  = zeros(Nz,1);
Nu_C  = zeros(Nz,1);
Cp_C  = zeros(Nz,1);

Re    = zeros(Nz,1);
Pr_LF = zeros(Nz,1);
Nu_LF = zeros(Nz,1);

h_L   = zeros(Nz,1);
h_V   = zeros(Nz,1);
h_G   = zeros(Nz,1);
h_G_s = zeros(Nz,1); %#ok<NASGU>
phi_T = zeros(Nz,1);
phi_T_x = ones(Nz,1);

q_r   = zeros(Nz,1);     % kW/m^2
q_CS  = zeros(Nz,1);
q_CL  = zeros(Nz,1);
q_LE  = zeros(Nz,1);

h_nb    = zeros(Nz,1);
h_conv  = zeros(Nz,1);
h_conv_1 = zeros(Nz,1);
h_conv_2 = zeros(Nz,1);
h_nb_2   = zeros(Nz,1);
h_conv_3 = zeros(Nz,1);
h_nb_3   = zeros(Nz,1);
Flw_4    = zeros(Nz,1);
Slw_4    = zeros(Nz,1);
h        = zeros(Nz,1);   % final HTC, kW/m^2-K

%% -----------------------------
% 2) Boundary conditions
% ------------------------------
G_V(1) = s.G_V0;

if G_V(1) == 0 || ~s.calculate_entraintment_0
    G_ED(1) = 1e-6;
else
    % Barbosa-Hewitt 2002 initial entrainment fraction form
    ratio = sqrt((rhoL*max(s.G_L0,1e-12)) / max(rhoV*G_V(1), 1e-12)) * Di^2;
    coeff = (1/(0.95e-2 + 342.55e-2*ratio) - 1);
    G_ED(1) = max(s.G_L0,1e-12) / max(coeff, 1e-12);
end

T_C(1) = T_BP(1);

% --- modified: G_L(1) is obtained from delta(1)
G_L(1) = calc_GL_from_delta( ...
    delta(1), G_V(1), G_ED(1), Di, ...
    rhoL, rhoV, rhoED, muL, muV, sigma, ...
    m.delta_calculation);

%% -----------------------------
% 3) Marching along z
% ------------------------------
for i = 1:Nz

    % ---- current local mass-based quantities ----
    u_L(i)  = G_L(i)  / rhoL;
    u_V(i)  = G_V(i)  / rhoV;
    u_ED(i) = G_ED(i) / rhoED;
    m_L(i) = G_L(i) * A;
    m_V(i) = G_V(i) * A;
    m_ED(i) = G_ED(i) * A;
    % Core averaged density
    rho_C(i) = (G_ED(i) + G_V(i)) / (G_ED(i)/rhoED + G_V(i)/rhoV + 1e-20);

    % Core velocity
    core_void_factor = max(1e-8, (1 - 2*delta(i)/Di)^2);
    u_C(i) = (G_ED(i)/rhoED + G_V(i)/rhoV) / core_void_factor;

    Re_C(i) = rho_C(i)*u_C(i)*Di / muV;
    Re_G(i) = rhoV*u_V(i)*Di / muV;
    We_G(i) = rhoV*Di*u_V(i)^2 / sigma;

    % Droplet concentration in gas core
    C_ED(i) = m_ED(i) / (m_ED(i)/rhoED + m_V(i)/rhoV + 1e-20);

    % Critical liquid flux for entrainment onset
    G_C(i) = (muL/Di) * exp(5.8504 + 0.4249*(muV/muL)*sqrt(rhoL/rhoV));

    % ---- deleted: all delta calculation content has been removed ----

    % Enforce geometric bound only
    delta(i) = max(1e-8, min(0.49*Di, delta(i)));

    % Recompute with prescribed delta
    Re_LF(i) = rhoL*u_L(i)*(4*delta(i)) / muL;
    f_SC(i)  = 0.079 / (abs(Re_C(i) + 1e-6)^0.25);
    f_I(i)   = f_SC(i) * (1 + 24*(delta(i)/Di)*(rhoL/rhoV)^(1/3));
    tau_I(i) = 0.5*rho_C(i)*f_I(i)*(u_C(i) - u_L(i))^2;
    f_LF(i)  = 0.079 / (abs(Re_LF(i) + 1e-6)^0.25);

    % ---- deposition coefficient k_D ----
    if C_ED(i)/rhoV <= 0.3
        k_D(i) = 0.18 * sqrt(sigma/(rhoV*Di));
    else
        k_D(i) = 0.083 * (rhoV/max(C_ED(i),1e-20))^0.65 * sqrt(sigma/(rhoV*Di));
    end

    % ---- entrainment source ED_cal ----
    if G_L(i) > G_C(i)
        ED_cal(i) = 5.75e-5 * G_V(i) * ...
            ((((rhoL*Di)/(sigma*rhoV^2))*((G_L(i)-G_C(i))^2))^0.316);
    else
        ED_cal(i) = 0.0;
    end

    % ---- selector for entrainment / deposition ----
    switch lower(m.entraintment_calculation)
        case 'no_entraintment'
            D(i)  = 0.0;
            ED(i) = 0.0;
        case 'cal_entraintment'
            D(i)  = k_D(i)*C_ED(i);
            ED(i) = ED_cal(i);
        otherwise
            error('Unknown entraintment_calculation option.');
    end

    % ---- thermal dimensionless groups ----
    X_vap(i) = min(0.9999, max(1e-4, G_V(i)/(G_V(i)+G_L(i)+G_ED(i)+1e-20)));
    X_tt(i)  = (muL/muV)^0.1 * (rhoV/rhoL)^0.5 * (((1-X_vap(i))/X_vap(i))^0.9);

    Re(i)    = rhoL*u_L(i)*Di / muL;
    Pr_LF(i) = muL*(cpL*1e3) / kL;
    Pr_C(i)  = muV*(cpV*1e3) / kV;

    % Single-phase HTC references
    h_L(i) = (kL/Di) * (0.023*(max(Re(i),1e-12)^0.8)*(Pr_LF(i)^0.4)) / 1e3; % kW/m^2-K
    h_V(i) = (kV/Di) * (0.023*(max(Re_C(i),1e-12)^0.8)*(Pr_C(i)^0.4)) / 1e3;

    % ---- core heat transfer ----
    Nu_C(i) = 0.023*(abs(Re_C(i)+1e-6)^0.8)*(Pr_C(i)^0.4);
    h_G(i)  = Nu_C(i)*kV/Di/1e3;           % kW/m^2-K

    % Equilibrium model
    T(i) = T_BP(i);

    switch lower(m.energy_balance_calculation)
        case 'equilibrium'
            if i > 1
                T_C(i) = T_C(i-1);
            end
            q_LE(i) = 0;
            q_CL(i) = 0;
            phi_T_x(i) = 1.0;
            phi_T(i)   = 0.0;

        case 'non_equilibrium'
            phi_T_x(i) = 1.0;
            Cp_C(i) = (G_V(i)+G_ED(i)) / ...
                      (G_V(i)/max(cpV,1e-20) + G_ED(i)/max(cpED,1e-20) + 1e-20);

            q_LE(i) = (ED(i) - D(i))*cpED*(T(i) - T_C(i));
            q_CL(i) = h_G(i)*(phi_T(i) + phi_T_x(i))*(T(i) - T_C(i));

            if i < Nz
                dTcdz = L*4*(h_G(i)*(T(i)-T_C(i)) + q_LE(i)) / ...
                        ((G_V(i)+G_ED(i))*Cp_C(i)*Di + 1e-20);
                T_C(i+1) = T_C(i) + dz*dTcdz;
            end

            val = h_G(i)*(T(i)-T_C(i));
            phi_T(i) = log(val/(phi_T_x(i)+1e-6) + 1.0);

        otherwise
            error('Unknown energy_balance_calculation option.');
    end

    % ---- initial guess of local evaporation from wall heat ----
    q_CS(i) = q_r(i) - q_CL(i) - q_LE(i);
    V_L(i)   = max(0, q_r(i) / max(DHvap,1e-12));
    V_ED(i)  = 0.0;
%     if i == 1
%         q_r(i) = max(0, h_G(i)*(T_w(i)-T(i)));
%     else
%         q_r(i) = max(0, q_r(i-1));
%     end
% 
%     V_L(i)   = max(0, q_r(i) / max(DHvap,1e-12));
%    
%     q_CS(i)  = V_L(i)*DHvap;

    % ---- boiling / convection correlations ----

    % Correlation 0
    h_nb(i)   = 55*(P/Pc)^0.12 * (-log10(P/Pc))^(-0.55) * Mw^(-0.5) * (max(q_r(i),1e-12)^0.67);
    Nu_LF(i)  = 0.0133*(max(Re_LF(i),1e-12)^0.69)*(Pr_LF(i)^0.4);
    h_conv(i) = Nu_LF(i)*kL/max(delta(i),1e-12)/1e3;

    % Correlation 1
    Fr_1 = ((G_L(i)/rhoL)^2)/(g*Di);
    if Fr_1 >= 0.25
        Fw_1 = 1.0;
    else
        Fw_1 = 1.32*(Fr_1^0.2);
    end
    h_conv_1(i) = Fw_1*h_L(i)*(1 + 1.925*(X_tt(i)^(-0.83)));

    % Correlation 2
    if X_tt(i) >= 9.979
        Fc_2 = 1.0;
    else
        Fc_2 = 2.35*(0.213 + 1/X_tt(i))^0.736;
    end
    h_conv_2(i) = Fc_2*h_L(i);
    Sc_2 = 1/(1 + 2.53e-6*(Re(i)*(Fc_2^1.25))^1.17);

    DT_sat_2 = 10;
    P_w = 10^(A_ant - B_ant/(C_ant + T(i) + DT_sat_2));
    DP_sat_2 = (P_w - P)*101325.0;

    h_nb_2(i) = 0.00122 * ...
        ((kL^0.79)*((1e3*cpL)^0.45)*(rhoL^0.49)) / ...
        ((sigma^0.5)*(muL^0.29)*(rhoV^0.24)*((1e3*DHvap)^0.24)) * ...
        (DT_sat_2^0.24) * (max(DP_sat_2,0)^0.75);

    % Correlation 3
    Fs_3 = ( ...
        ((1-X_vap(i))^1.5 + 1.9*(X_vap(i)^0.6)*((1-X_vap(i))^0.01)*(rhoL/rhoV)^0.35)^(-2.2) ...
        + (h_V(i)/h_L(i))*(X_vap(i)^0.01)*(1 + 8*((1-X_vap(i))^0.7)*(rhoL/rhoV)^0.67)^(-2.0) ...
        )^(-0.5);

    h_conv_3(i) = Fs_3*h_L(i);

    Cf_3  = 0.435*((Mw/2.016)^0.27);
    nf_3  = 0.8 - 0.1*(10^(0.76*(P/Pc)));
    dPsatdT = log10(exp(1))*P*B_ant / ((C_ant + T(i))^2);
    q_crit_3 = 2*sigma*T(i)*h_L(i)/(3*rhoV*DHvap);
    h_pb = 3.58*(dPsatdT/sigma)^0.6;

    if q_r(i) <= q_crit_3
        h_nb_3(i) = 0;
    else
        h_nb_3(i) = h_pb*Cf_3*(q_r(i)/20)^nf_3;
    end

    % Correlation 4
    denom = rhoV - 1.0;
    if (1 + X_vap(i)*Pr_LF(i)*(rhoL/denom)) <= 0
        Flw_4(i) = (1 + X_vap(i)*Pr_LF(i)*(rhoL/(rhoV + 0.0)))^0.35;
    else
        Flw_4(i) = (1 + X_vap(i)*Pr_LF(i)*(rhoL/denom))^0.35;
    end
    Slw_4(i) = 1/(1 + 0.055*(Flw_4(i)^0.1)*(Re(i)^0.16));

    % ---- final correlation selector ----
    switch lower(m.h_calculation)
        case 'simple'
            h(i) = 300.0/1e3;

        case 'correlation_0'
            h(i) = (h_nb(i)^3 + h_conv(i)^3)^(1/3);

        case 'correlation_1'
            h(i) = (h_nb(i)^2.5 + h_conv_1(i)^2.5)^(1/2.5);

        case 'correlation_2'
            h(i) = h_conv_2(i) + Sc_2*h_nb_2(i);

        case 'correlation_3'
            h(i) = (h_conv_3(i)^3 + h_nb_3(i)^3)^(1/3);

        case 'correlation_4'
            h(i) = ((Slw_4(i)*h_nb(i))^2 + (Flw_4(i)*h_L(i))^2)^(1/2);

        case 'correlation_a'
            h(i) = ( ...
                (h_nb(i)^2.5 + h_conv_1(i)^2.5)^(1/2.5) + ...
                h_conv_2(i) + Sc_2*h_nb_2(i) + ...
                (h_conv_3(i)^3 + h_nb_3(i)^3)^(1/3) + ...
                ((Slw_4(i)*h_nb(i))^2 + (Flw_4(i)*h_L(i))^2)^(1/2) ) / 4.0;

        otherwise
            error('Unknown h_calculation option.');
    end

    % ---- wall heat flux ----
    q_r(i) = h(i) * (T_w(i) - T(i));
    q_r(i) = max(q_r(i), 0);

    % update evaporation term with corrected q_r
    V_L(i)  = max(0, q_r(i) / max(DHvap,1e-12));
    q_CS(i) = V_L(i)*DHvap;

    % ---- march mass balances to next cell ----
    if i < Nz
        dG_L_dz  = (L*4*(D(i) - ED(i) - V_L(i))) / Di;
        dG_V_dz  = (L*4*(V_L(i) + V_ED(i))) / Di;
        dG_ED_dz = (L*4*(ED(i) - D(i) - V_ED(i))) / Di;

        if G_L(i) == 0
            dG_L_dz = (L*4*max(0, D(i) - ED(i) - V_L(i))) / Di;
        end
        if G_ED(i) == 0
            dG_ED_dz = (L*4*max(0, ED(i) - D(i) - V_ED(i))) / Di;
        end

        G_V(i+1)  = max(0, G_V(i)  + dz*dG_V_dz);
        G_ED(i+1) = max(0, G_ED(i) + dz*dG_ED_dz);

        % modified: G_L(i+1) obtained from prescribed delta(i+1)
        G_L(i+1) = calc_GL_from_delta( ...
            delta(i+1), G_V(i+1), G_ED(i+1), Di, ...
            rhoL, rhoV, rhoED, muL, muV, sigma, ...
            m.delta_calculation);
    end
end

%% -----------------------------
% 4) Total heat transfer
% ------------------------------
dA = pi*Di*(L*dz);
Q_total = sum(q_r*1e3 * dA);   % [W], q_r was kW/m^2

%% -----------------------------
% 5) Output
% ------------------------------
out.z     = z*L;
out.q     = q_r;        % kW/m^2
out.Q     = Q_total;    % W
out.h     = h;          % kW/m^2-K
out.delta = delta;
out.G_L   = G_L;
out.G_V   = G_V;
out.G_ED  = G_ED;
out.T     = T;
out.T_C   = T_C;
out.V_L   = V_L;
out.D     = D;
out.ED    = ED;
out.tau_I = tau_I;

end


%% ============================================================
function G_L = calc_GL_from_delta(delta, G_V, G_ED, Di, ...
    rhoL, rhoV, rhoED, muL, muV, sigma, delta_mode)
%CALC_GL_FROM_DELTA
% Output: G_L
%
% 1) triangular_rel:
%    (G_L^2)*(f_LF) = ((4*delta/Di)^2)*(2*tau_I*rhoL)
%
% 2) vol_fraction:
%    delta = 0.5*Di*w_frac_L
%    where
%    w_frac_L = (G_L/rhoL) / (G_L/rhoL + G_V/rhoV + G_ED/rhoED)

delta = max(1e-8, min(0.49*Di, delta));

switch lower(delta_mode)

    case 'triangular_rel'
        G_L = solve_GL_triangular(delta, G_V, G_ED, Di, ...
            rhoL, rhoV, rhoED, muL, muV, sigma);

    case 'vol_fraction'
        beta = 2*delta/Di;   % w_frac_L
        beta = max(1e-12, min(0.999999, beta));

        S = G_V/rhoV + G_ED/rhoED;
        % beta = x / (x + S), x = G_L/rhoL
        x = beta*S / max(1-beta, 1e-12);
        G_L = rhoL * x;

    otherwise
        error('Unknown delta_calculation option.');
end

G_L = max(0, G_L);

end


%% ============================================================
function G_L = solve_GL_triangular(delta, G_V, G_ED, Di, ...
    rhoL, rhoV, rhoED, muL, muV, sigma)
% Solve for G_L from:
% (G_L^2)*(f_LF) = ((4*delta/Di)^2)*(2*tau_I*rhoL)

rhoC = (G_ED + G_V) / (G_ED/rhoED + G_V/rhoV + 1e-20);
uC   = (G_ED/rhoED + G_V/rhoV) / max(1e-8,(1 - 2*delta/Di)^2);
ReC  = rhoC*uC*Di / muV;
fSC  = 0.079/(abs(ReC + 1e-6)^0.25);

% initial guess
G_L = max(1e-6, rhoL*sqrt(max(1e-12, sigma/(rhoV*Di))));

for it = 1:80
    uL   = G_L / rhoL;
    ReLF = rhoL*uL*(4*delta)/muL;
    fLF  = 0.079/(abs(ReLF + 1e-6)^0.25);

    fI   = fSC*(1 + 24*(delta/Di)*(rhoL/rhoV)^(1/3));
    tauI = 0.5*rhoC*fI*(uC - uL)^2;

    F = (G_L^2)*fLF - ((4*delta/Di)^2)*(2*tauI*rhoL);

    dG = max(1e-8, 1e-5*max(G_L,1));
    Gp = G_L + dG;

    uLp   = Gp / rhoL;
    ReLFp = rhoL*uLp*(4*delta)/muL;
    fLFp  = 0.079/(abs(ReLFp + 1e-6)^0.25);
    tauIp = 0.5*rhoC*fI*(uC - uLp)^2;

    Fp = (Gp^2)*fLFp - ((4*delta/Di)^2)*(2*tauIp*rhoL);
    dF = (Fp - F)/dG;

    G_new = G_L - F/(dF + 1e-20);
    G_new = max(1e-8, G_new);

    if abs(G_new - G_L) < 1e-8
        G_L = G_new;
        return;
    end

    G_L = G_new;
end

end