% [x] = hw8solve(A, b)
%
% Return a global minimizer for the constrained optimization problem
%   minimize x'Ax/2 - x'b s.t. x'x = 1
% Your code should run in an overall time of O(n^3).
%
function [x] = hw8solve(A, b)