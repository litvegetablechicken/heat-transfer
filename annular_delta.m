function out = annular_delta(param, q_r,Boundary,delta,L)
% ANNULAR_HEAT_TRANSFER
% Grid-based MATLAB implementation of annular-region equations.
%
% INPUT
%   param : struct
%   q_r   : wall heat flux [kW/m^2], Nz x 1
%
% OUTPUT
%   out   : struct

%% -------------------------------------------------
% 0) Basic checks
%% -------------------------------------------------

Nz = param.geom.NzA;
Di = param.geom.Di;

if isrow(q_r)
    q_r = q_r(:);
end

assert(length(q_r)==Nz,'q_r must have length Nz')

z  = linspace(0,1,Nz)';
dz = z(2)-z(1);

A = pi*(Di/2)^2;

f = param.fluid;
model = param.model;

rhoL = f.rhoL;
rhoV = f.rhoV;
rhoED = f.rhoED;

muL = f.muL;
muV = f.muV;

kL = f.kL;
kV = f.kV;

cpL = f.cpL;
cpV = f.cpV;

sigma = f.sigma;
DHvap = f.DH_vap;

P  = f.P;
Pc = f.Pc;
Mw = f.Mw;

A_ant = f.A_ant;
B_ant = f.B_ant;
C_ant = f.C_ant;

if isscalar(f.Tbp)
    T_BP = f.Tbp*ones(Nz,1);
else
    T_BP = f.Tbp(:);
end

g = 9.81;

%% -------------------------------------------------
% 1) Allocate variables
%% -------------------------------------------------
G_L  = zeros(Nz,1);
G_V  = zeros(Nz,1);
G_ED = zeros(Nz,1);

T     = T_BP;
T_C   = zeros(Nz,1);
T_w   = zeros(Nz,1);


u_L  = zeros(Nz,1);
u_V  = zeros(Nz,1);
u_ED = zeros(Nz,1);

m_L = zeros(Nz,1);
m_V = zeros(Nz,1);
m_ED = zeros(Nz,1);

rho_C = zeros(Nz,1);
u_C   = zeros(Nz,1);

Re_C = zeros(Nz,1);
Re_G = zeros(Nz,1);
We_G = zeros(Nz,1);

tau_I = zeros(Nz,1);
f_SC  = zeros(Nz,1);
f_I   = zeros(Nz,1);
f_LF  = zeros(Nz,1);
Re_LF = zeros(Nz,1);

C_ED  = zeros(Nz,1);
k_D   = zeros(Nz,1);
G_C   = zeros(Nz,1);
ED_cal = zeros(Nz,1);

D  = zeros(Nz,1);
ED = zeros(Nz,1);

V_L  = zeros(Nz,1);
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

phi_T = zeros(Nz,1);
phi_T_x = ones(Nz,1);

q_CS = zeros(Nz,1);
q_CL = zeros(Nz,1);
q_LE = zeros(Nz,1);

h_nb    = zeros(Nz,1);
h_conv  = zeros(Nz,1);

h_conv_1 = zeros(Nz,1);
h_conv_2 = zeros(Nz,1);
h_nb_2   = zeros(Nz,1);

h_conv_3 = zeros(Nz,1);
h_nb_3   = zeros(Nz,1);
G_LF = zeros(Nz,1);
Flw_4 = zeros(Nz,1);
Slw_4 = zeros(Nz,1);
F = zeros(Nz,1);
h = zeros(Nz,1);

%% -------------------------------------------------
% 2) Boundary conditions
%% -------------------------------------------------

G_V(1) = Boundary.A.G_V;
G_L(1)=Boundary.A.G_L;
if G_V(1)==0 || ~model.calculate_entraintment_0
    G_ED(1)=1e-6;
else

    ratio = sqrt((rhoL*max(G_L(1),1e-12))/(rhoV*G_V(1)))*Di^2;

    coeff = (1/(0.95e-2+342.55e-2*ratio)-1);

    G_ED(1)=max(G_L(1),1e-12)/max(coeff,1e-12);

end

T_C(1)=T_BP(1);



%% -------------------------------------------------
% 3) Marching along z
%% -------------------------------------------------

for i=1:Nz

    u_L(i)=G_L(i)/rhoL;
    u_V(i)=G_V(i)/rhoV;
    u_ED(i)=G_ED(i)/rhoED;

    m_L(i)=G_L(i)*A;
    m_V(i)=G_V(i)*A;
    m_ED(i)=G_ED(i)*A;

    rho_C(i)=(G_ED(i)+G_V(i))/(G_ED(i)/rhoED+G_V(i)/rhoV+1e-20);
    w_frac_L = m_L(i) / rhoL / (m_V(i)/rhoV + m_ED(i) / rhoED +  m_L(i) / rhoL );
    % 计算通过三角关系计算液膜厚度
    core_void_factor=max(1e-8,(1-2*delta(i)/Di)^2);
    u_C(i)=(G_ED(i)/rhoED+G_V(i)/rhoV)/core_void_factor;
    Re_C(i)=rho_C(i)*u_C(i)*Di/muV;
    Re_G(i)=rhoV*u_V(i)*Di/muV;
    We_G(i)=rhoV*Di*u_V(i)^2/sigma;

    %% triangular relationship
    Re_LF(i)=rhoL*u_L(i)*(4*delta(i))/muL;

    f_SC(i)=0.079/(abs(Re_C(i)+1e-6)^0.25);

    f_I(i)=f_SC(i)*(1+24*(delta(i)/Di)*(rhoL/rhoV)^(1/3));

    tau_I(i)=0.5*rho_C(i)*f_I(i)*(u_C(i)-u_L(i))^2;
    if Re_LF(i) > 1e3
    f_LF(i)=0.079/(abs(Re_LF(i)+1e-6)^0.25);
    else
            f_LF(i)=16/Re_LF(i);

    end
    G_LF(i) = 4/Di * delta(i) * sqrt(2 * tau_I(i) * rhoL / f_LF(i));
    %% deposition and entrainment fluxes
    C_ED(i)=m_ED(i)/(m_ED(i)/rhoED+m_V(i)/rhoV+1e-20);
    if C_ED(i)/rhoV<=0.3
        k_D(i)=0.18*sqrt(sigma/(rhoV*Di));
    else
        k_D(i)=0.083*(rhoV/max(C_ED(i),1e-20))^0.65*sqrt(sigma/(rhoV*Di));
    end
    G_C(i)=(muL/Di)*exp(5.8504+0.4249*(muV/muL)*sqrt(rhoL/rhoV));

    if G_L(i)>G_C(i)
        ED_cal(i)=5.75e-5*G_V(i)*...
            ((((rhoL*Di)/(sigma*rhoV^2))*((G_L(i)-G_C(i))^2))^0.316);
    else
        ED_cal(i)=0;
    end
    if model.calculate_entraintment_0 == true
        D(i)=k_D(i)*C_ED(i);
        ED(i)=ED_cal(i);
    else
        D(i)=0;
        ED(i)=0;
    end

    %% heat transfer correlations
    Re(i)=rhoL*u_L(i)*Di/muL;
    Pr_LF(i)=muL*(cpL)/kL;
    Pr_C(i)=muV*(cpV)/kV;

 

    Nu_C(i)=0.023*(abs(Re_C(i)+1e-6)^0.8)*(Pr_C(i)^0.4);

    h_G(i)=Nu_C(i)*kV/Di/1e3;

    T(i)=T_BP(i);
    
    X_vap(i)=min(0.9999,max(1e-4,G_V(i)/(G_V(i)+G_L(i)+G_ED(i)+1e-20)));
    X_tt(i)=(muL/muV)^0.1*(rhoV/rhoL)^0.5*(((1-X_vap(i))/X_vap(i))^0.9);
    switch lower(model.h_calculation)
        case 'correlation_0'
            % Correlation 0 Chen [16]
            h_nb(i)=55*(P/Pc)^0.12*(-log10(P/Pc))^(-0.55)*Mw^(-0.5)*(max(q_r(i),1e-12)^0.67);

            Nu_LF(i)=0.0133*(max(Re_LF(i),1e-12)^0.69)*(Pr_LF(i)^0.4);
            h_conv(i)=Nu_LF(i)*kL/max(delta(i),1e-12)/1e3;
            h_nb(i) = h_nb(i) /1e3;
            h(i)=((h_nb(i)^3+h_conv(i)^3)^(1/3))*1e3;
        case 'correlation_1'
            % Correlation 1 Cooper [17]
            h_L(i)=(kL/Di)*(0.023*(max(Re(i),1e-12)^0.8)*(Pr_LF(i)^0.4));
            h_V(i)=(kV/Di)*(0.023*(max(Re_C(i),1e-12)^0.8)*(Pr_C(i)^0.4));
            Fr_1 = ((G_L(i)/rhoL)^2)/(g*Di);
            if Fr_1 >= 0.25
                Fw_1 = 1.0;
            else
                Fw_1 = 1.32*(Fr_1^0.2);
            end
            h_nb(i) = h_nb(i) /1e3;
            h_conv_1(i) = Fw_1*h_L(i)*(1 + 1.925*(X_tt(i)^(-0.83)));
            h_conv_1(i) = h_conv_1(i) /1e3;
            h(i) = ((h_nb(i)^2.5 + h_conv_1(i)^2.5)^(1/2.5))*1e3;
        case 'correlation_2' % Chen approximation [18]
            if X_tt(i) >= 9.979
                Fc_2 = 1.0;
            else
                Fc_2 = 2.35*(0.213 + 1/X_tt(i))^0.736;
            end
            h_L(i)=(kL/Di)*(0.023*(max(Re(i),1e-12)^0.8)*(Pr_LF(i)^0.4));

            h_conv_2(i) = Fc_2*h_L(i);
            Sc_2 = 1/(1 + 2.53e-6*(max(Re(i),1e-12)*(Fc_2^1.25))^1.17);
            fun = @(x)Tw_balance(x,T(i),h_conv_2(i),Sc_2,f,q_r(i));
            options = optimoptions('fsolve', 'Display','off');
            x0 = T(i);
            [Tw,~,exitflag] = fsolve(fun,x0,options);

            if exitflag <0
                % error('error of delta calculation')
                disp('error of correlation_2')
            elseif exitflag == 0
                disp('MaxIterations exceeds MaxFunctionEvaluations')
            end

            DT_sat_2 = Tw - T(i);

            P_w = refpropm('P','T', Tw,'Q',0,f.name)*1e2;% 输出kPa

            DP_sat_2 = (P_w - P*1e5);
            h_nb_2(i) = 0.00122 * ...
                ((kL^0.79)*((cpL)^0.45)*(rhoL^0.49)) / ...
                ((sigma^0.5)*(muL^0.29)*(rhoV^0.24)*((DHvap)^0.24)) * ...
                (DT_sat_2^0.24) * (max(DP_sat_2,0)^0.75);
            h(i) = h_conv_2(i) + Sc_2*h_nb_2(i);
        case 'correlation_3' % Steiner [20]
            h_L(i)=(kL/Di)*(0.023*(max(Re(i),1e-12)^0.8)*(Pr_LF(i)^0.4));
            h_V(i)=(kV/Di)*(0.023*(max(Re_C(i),1e-12)^0.8)*(Pr_C(i)^0.4));
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
                h_nb_3(i) = h_pb*Cf_3*(q_r(i)/1e3/20)^nf_3;
            end
            h_conv_3(i) = h_conv_3(i)/1e3;
            h_nb_3(i) = h_nb_3(i);
            h(i) = (h_conv_3(i)^3 + h_nb_3(i)^3)^(1/3);
            h(i) = h(i) * 1e3;
        case 'correlation_4'
            denom = rhoV - 1.0;
            h_nb(i)=55*(P/Pc)^0.12*(-log10(P/Pc))^(-0.55)*Mw^(-0.5)*(max(q_r(i),1e-12)^0.67);
            h_nb(i) = h_nb(i) / 1e3;
            h_L(i)=(kL/Di)*(0.023*(max(Re(i),1e-12)^0.8)*(Pr_LF(i)^0.4));
            h_L(i) = h_L(i)/1e3;
            if (1 + X_vap(i)*Pr_LF(i)*(rhoL/denom)) <= 0
                Flw_4(i) = (1 + X_vap(i)*Pr_LF(i)*(rhoL/(rhoV + 0.0)))^0.35;
            else
                Flw_4(i) = (1 + X_vap(i)*Pr_LF(i)*(rhoL/denom))^0.35;
            end
            Slw_4(i) = 1/(1 + 0.055*(Flw_4(i)^0.1)*(Re(i)^0.16));
            h(i) = ((Slw_4(i)*h_nb(i))^2 + (Flw_4(i)*h_L(i))^2)^(1/2);
            h(i) = h(i) * 1e3;

    end
    
    %% wall temperature from heat flux
    T_w(i)=T(i)+q_r(i)/max(h(i),1e-12);
    %% evaporation
    V_L(i)=max(0,q_r(i)/max(DHvap,1e-12));

    q_CS(i)=V_L(i)*DHvap;

    %% -------------------------------------------------
    % march to next cell
    %% -------------------------------------------------

    if i<Nz

        dG_L_dz=(L*4*(D(i)-ED(i)-V_L(i)))/Di;

        dG_V_dz=(L*4*(V_L(i)+V_ED(i)))/Di;

        dG_ED_dz=(L*4*(ED(i)-D(i)-V_ED(i)))/Di;

        G_V(i+1)=G_V(i)+dz*dG_V_dz;

        G_ED(i+1)=G_ED(i)+dz*dG_ED_dz;

        G_L(i+1)=G_L(i)+dz*dG_L_dz;
        F(i) =  cal_delta(delta(i),Di,G_ED(i),rhoED,G_V(i),rhoV,rhoL,u_L(i),muL,muV,rho_C(i),G_L(i));
    end

end

%% -------------------------------------------------
% total heat transfer
%% -------------------------------------------------

dA=pi*Di*(L*dz);

Q_total=sum(q_r*dA);

%% -------------------------------------------------
% output
%% -------------------------------------------------
out.term = G_L(end);
out.z=z*L;
out.V_L = V_L;
out.q=q_r;
out.D = D;
out.ED = ED;
out.Q=Q_total;

out.h=h;

out.delta=delta;

out.G_L=G_L;

out.G_V=G_V;

out.G_ED=G_ED;

out.T=T;

out.T_C=T_C;

out.Tw=T_w;
out.F =F;
end


%% ==========================================================
function F = cal_delta(x,Di,G_ED,rhoED,G_V,rhoV,rhoL,u_L,muL,muV,rho_C,G_L)
delta = x;
core_void_factor=max(1e-8,(1-2*delta/Di)^2);

u_C=(G_ED/rhoED+G_V/rhoV)/core_void_factor;
Re_LF=rhoL*u_L*(4*delta)/muL;
Re_C=rho_C*u_C*Di/muV;
f_SC=0.079/(abs(Re_C+1e-6)^0.25);

f_I=f_SC*(1+24*(delta/Di)*(rhoL/rhoV)^(1/3));

tau_I=0.5*rho_C*f_I*(u_C-u_L)^2;
% if Re_LF > 1e3
%     f_LF=0.079/(abs(Re_LF+1e-6)^0.25);
% else
%     f_LF=16/Re_LF;
% 
% end
    f_LF=0.079/(abs(Re_LF+1e-6)^0.25);
G_LF = 4/Di * delta * sqrt(2 * tau_I * rhoL / f_LF);
F = G_L - G_LF;
end


function F = Tw_balance(x,T,h_conv_2,Sc_2,f,q_r)
Tw = x;
DT_sat_2 = Tw - T;

P_w = refpropm('P','T', Tw,'Q',0,f.name)/1e2; % bar

DP_sat_2 = (P_w - f.P)*1e5;
h_nb_2 = 0.00122 * ...
    ((f.kL^0.79)*((f.cpL)^0.45)*(f.rhoL^0.49)) / ...
    ((f.sigma^0.5)*(f.muL^0.29)*(f.rhoV^0.24)*((f.DHvap)^0.24)) * ...
    (DT_sat_2^0.24) * (max(DP_sat_2,0)^0.75);
h = h_conv_2 + Sc_2*h_nb_2;
T_w=T+q_r/max(h,1e-12);
F = Tw - T_w;
end