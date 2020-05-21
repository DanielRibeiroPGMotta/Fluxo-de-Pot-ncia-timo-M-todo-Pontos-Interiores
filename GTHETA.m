% GTHETA computes the Hessian matrix of the
% line flow functions and modifies it to account for the
% reduction of the coefficient matrix
%
function Gt = GTHETA(na,nb,b,theta,yl_mag,ang_y,A2,lambda,wubar,wbar,F2,rubar,rbar)
%
% Computing zeta
%
zeta = A2' * lambda + wubar + wbar;
%
% Computing Gtheta, 1st step
%
[n,m] = size(theta);
[nl,m] = size(na);
n = n - 1;
Gt = zeros(n,n);
cos_t_a1 = cos(theta(na)-theta(nb)+ang_y);
cos_t_a2 = cos(theta(nb)-theta(na)+ang_y);
bz1 = yl_mag .* zeta(1:nl);
bz2 = yl_mag .* zeta(nl+1:2*nl);
bzs1 = bz1 .* cos_t_a1;
bzs2 = bz2 .* cos_t_a2;
%
for k = 1 : nl
    na1 = na(k) - 1;
    nb1 = nb(k) - 1;
    if na(k) > 1
       Gt(na1,na1) = Gt(na1,na1) - bzs1(k) - bzs2(k);
       Gt(nb1,na1) = Gt(nb1,na1) + bzs1(k) + bzs2(k);
       Gt(na1,nb1) = Gt(na1,nb1) + bzs1(k) + bzs2(k);
    end
    Gt(nb1,nb1) = Gt(nb1,nb1) - bzs1(k) - bzs2(k);
end
%
% Computing Gtheta, 2nd step
%
M = rubar .\ wubar - rbar .\ wbar;
Gt = Gt - F2' * diag(M) * F2;
