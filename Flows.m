function [f,F2]=flows(n,nl,na,nb,b,theta,yl_mag,ang_y)
thetaij = theta(na)-theta(nb);
%bdiag = diag(b);
%f = bdiag * sin(thetaij);
f(1:nl,1) = yl_mag .* (cos(thetaij+ang_y) - cos(ang_y));
f(nl+1:2*nl,1) = yl_mag .* (cos(-thetaij+ang_y) - cos(ang_y));
%
F2 = zeros(2*nl,n-1);
for i = 1: nl
%    costheta = cos (thetaij(i));
	 sin_t_a1 = sin(thetaij(i)+ang_y(i));
         sin_t_a2 = sin(-thetaij(i)+ang_y(i));
    if na(i) > 1, 
       %F(i,na(i)-1) = costheta; 
        F2(i,na(i)-1) = -yl_mag(i)* sin_t_a1; 
        F2(i+nl,na(i)-1) = yl_mag(i)* sin_t_a2;
    end 
    %F(i,nb(i)-1) = - costheta; 
	 F2(i,nb(i)-1) = yl_mag(i)*sin_t_a1;
         F2(i+nl,nb(i)-1) = -yl_mag(i)*sin_t_a2;
end
%F = bdiag * F2;
%F = diag(yl_mag) * F2;
