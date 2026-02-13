
function [lzl,lzsa,lza,lzm,lzv,Torc,Tex,h,Tw_d,Tw_do] = datareduction(name)
excel = name;
[lzl,lzsa,lza,lzm,lzv] = readL(excel,'HEX_DP_5FP.Liquid.L','HEX_DP_5FP.Annular.L','HEX_DP_5FP.Vapor.L','HEX_DP_5FP.Subannular.L','HEX_DP_5FP.Mist.L');

Torc = readata(excel,'HEX_DP_5FP.Liquid.T','HEX_DP_5FP.Annular.T','HEX_DP_5FP.Vapor.T',lzl,lza,lzv,'HEX_DP_5FP.Subannular.T','HEX_DP_5FP.Mist.T',lzsa,lzm);
Torc(:,2) = Torc(:,2)-273.15;

Tex = readata(excel,'HEX_DP_5FP.Annulus_L.T_ex','HEX_DP_5FP.Annulus_A.T_ex','HEX_DP_5FP.Annulus_V.T_ex',lzl,lza,lzv,'HEX_DP_5FP.Annulus_SA.T_ex','HEX_DP_5FP.Annulus_M.T_ex',lzsa,lzm);
Tex(:,2) = Tex(:,2)-273.15;
Tw_d = readata(excel,'HEX_DP_5FP.Wall_L.Tw_d','HEX_DP_5FP.Wall_A.Tw_d','HEX_DP_5FP.Wall_V.Tw_d',lzl,lza,lzv,'HEX_DP_5FP.Wall_SA.Tw_d','HEX_DP_5FP.Wall_M.Tw_d',lzsa,lzm);
Tw_do = readata(excel,'HEX_DP_5FP.Wall_L.Tw_do','HEX_DP_5FP.Wall_A.Tw_do','HEX_DP_5FP.Wall_V.Tw_do',lzl,lza,lzv,'HEX_DP_5FP.Wall_SA.Tw_do','HEX_DP_5FP.Wall_M.Tw_do',lzsa,lzm);

h = readata(excel,'HEX_DP_5FP.Liquid.h','HEX_DP_5FP.Annular.h','HEX_DP_5FP.Vapor.h',lzl,lza,lzv,'HEX_DP_5FP.Subannular.h','HEX_DP_5FP.Mist.h',lzsa,lzm);
end


% Twex = [0.03545	31.74332
% 0.48157	31.99289
% 0.95018	32.7294
% 1.41875	33.34425
% 1.87635	34.44558
% 2.334	35.66856
% 2.81404	37.13518
% 3.26082	39.08796
% 3.73023	41.89266
% 4.18872	45.3055
% 4.64744	49.32663
% 5.11756	53.95621
% 5.5752	55.1792];
% x = Tex(:,1);
% y = Tex(:,2);
% xq = Twex(:,1);
% yq = interp1(x,y,xq,'linear');
% error = abs(yq - Twex(:,2))./Twex(:,2);

