% Return an approximate eigenvector v s.t. r = A*v-mu*v is minimal (in norm)
% where v is from the space spanned by V (V'*V = I).  Also return the
% norm of the min backward error E s.t. (A+E)*v = mu*v.
%
function [v, normE] = p3eigen(A, mu, V)