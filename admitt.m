% ADMITT computes the magnitudes and phases
% of off-diagonal elements of YBUS matrix
%
% Requires arrays RL and XL with the resistances
% and reactances of the transmission lines
%
icmplx = sqrt(-1);
den = RL .^ 2 + XL .^ 2;
yl = (RL - icmplx * XL) ./ den;
yl_mag = abs(yl);
ang_y = - (angle(yl) + pi);

