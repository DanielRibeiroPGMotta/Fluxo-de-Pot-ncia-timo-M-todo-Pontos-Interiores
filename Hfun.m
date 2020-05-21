% HFUN computes the (symmetrical of the) third and last segment of
% the rhs for the reduced linear system
%
function h = HFUN(A,Ag,p,f,d)
%
h = - Ag * p + A * f + d;
