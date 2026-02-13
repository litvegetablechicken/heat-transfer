function [Wt, Wp, Qeva, Qcon, H2, P2, H4, P4, H1] = ORC(fluid, T1, P1, Tcon, Tsc, Ets, m)
%本函数可以计算ORC的输出功率以及效率
%  working condition
% Thso K heat source temperature
% Tcso K cold source temperature
% Tap K approch temperature
% Ets isentropy efficiency
% m mass flow
% ORC parameters
% output
% Tsh turbine inlet superheat temperature
% Tsc condensor outlet supercooling temperature
Teva = refpropm('T','P',P1,'Q',0,fluid);
Tsh = T1 - Teva;
if Tsh == 0
    S1 = refpropm('S', 'T', T1, 'Q', 1,  fluid); % J/(kg K)
else
    S1 = refpropm('S', 'T', T1, 'P', P1,  fluid); % J/(kg K)
end
H1 = refpropm('H', 'P', P1, 'S', S1,  fluid); % J/kg
P2 = refpropm('P', 'T', Tcon, 'Q', 0, fluid); % kPa condensation temperature
H2s = refpropm('H', 'P', P2, 'S', S1, fluid);
H2 = H1 - (H1 - H2s) * Ets;
T2 = refpropm('T',  'P', P2, 'H', H2, fluid);

Wt = (H1 - H2) * m ;
% pump
P3 = P2;
if Tsc == 0
    H3 = refpropm('H', 'P', P3, 'Q', 0, fluid);
    S3 = refpropm('S', 'P', P3, 'Q', 0, fluid);
else
    T3 = Tcon - Tsc;
    H3 = refpropm('H', 'T', T3, 'P', P3, fluid);
    S3 = refpropm('S', 'T', T3, 'P', P3, fluid);
end
P4 = P1;
H4s = refpropm('H', 'P', P4, 'S', S3, fluid);
H4 = (H4s - H3) / Ets + H3;
Wp = (H4 - H3) * m;
Qeva = (H1 - H4) * m; % Qh heat absorb

Qcon = (H2 - H3) * m;  % Qc heat release

end