% TFUN computes the 1st segment of the rhs vector for the
% reduced linear system
%
function [t,t1]=TFUN(c,Q,p,Ag,lambda,piubar,pibar,subar,sbar,vubar,vbar)
%
gradC = c + Q * p;
t1 = gradC - Ag' * lambda + piubar + pibar;
%
t = - t1 + subar .\ vubar - sbar .\ vbar;
