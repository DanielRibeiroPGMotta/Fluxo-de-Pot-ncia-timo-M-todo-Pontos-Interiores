% UFUN computes the 2nd segment of the rhs vector for the
% reduced linear system
%
function [u,u1] = UFUN(F,A,lambda,wubar,wbar,rubar,rbar,yubar,ybar,zubar,zbar)
%
u1 = F' * (A' * lambda + wubar + wbar);
ywzubar = yubar + wubar .* zubar;
ywzbar = ybar + wbar .* zbar;
u =  - u1 + F' * (rubar .\ ywzubar - rbar .\ ywzbar);
