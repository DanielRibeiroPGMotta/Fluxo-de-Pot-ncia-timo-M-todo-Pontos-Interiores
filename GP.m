% GP computes the Hessian matrix of the cost function and
% modifies it to account for the reduction of the cofficient matrix
%
function G = GP(Q,subar,sbar,piubar,pibar)
%
G = Q - diag(subar.\piubar) + diag(sbar.\pibar);
