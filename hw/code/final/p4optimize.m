% Minimize x'*A*X s.t. x'*M*x = 1 and C*x = b.
% Return feasible=0 if no points satisfy the constraints.
%
function [x, feasible] = p4optimize(A, M, C, b)

