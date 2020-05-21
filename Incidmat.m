
% Bus-to-generator incidence matrix
Ag = zeros(n,ng);
for i = 1:ng
    Ag(ip(i),i) = 1;
end
% Bus-to-branch incidence matrix
A2 = zeros(n,nl);
for i = 1:nl
    A2(na(i),i) = 1;
    A2(nb(i),nl+i) = 1;
end
