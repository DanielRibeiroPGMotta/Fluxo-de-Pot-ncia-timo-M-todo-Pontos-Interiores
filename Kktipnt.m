% KKTIPNT checks the KKT conditions for the Interior Point
% Newton solution of the Economic Dispatch problem
%
DpL = c + Q * p -Ag' * lambda + piubar + pibar;
DthetaL =  F' * (  A2' * lambda + wubar + wbar);
DlambdaL = -Ag * p + A2 * f + d;
DpiubarL = p - subar - pubar;
DpibarL = p + sbar - pbar;
Dwubar = f - rubar -fubar;
Dwbar = f + rbar - fbar;
kkt1=[DpL;DthetaL;DlambdaL;DpiubarL;DpibarL;Dwubar;Dwbar];

