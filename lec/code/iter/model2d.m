% [T, u, f, h] = model2d(n)
%
% Set up
%    h^-2 * T * u = f
% where h^{-2} T is the discretized 2D Poisson problem on [0,1]^2
% with Dirichlet BCs using an (n+2)-by-(n+2) grid (including the
% constrained boundary nodes).
%
function [h, T, f, u] = model2d(n)

  % Construct the (sparse) 1D model matrix
  B = -ones(n,3);
  B(:,2) = 2;
  T1d = spdiags(B,[-1, 0, 1], n, n);

  % Construct the 2D version via Kronecker products
  I1d = speye(n);
  T = kron(I1d, T1d) + kron(T1d, I1d);
  h = 1/(n+1);

  % Build interior mesh (1D domain)
  x = linspace(0,1,n+2).';
  x = x(2:end-1);

  % Build 2D interior mesh
  [xx, yy] = meshgrid(x, x);
  xx = reshape(xx, n*n, 1);
  yy = reshape(yy, n*n, 1);

  % Constant right hand side + true solution via sparse direct
  if nargout > 2, f = ones(n*n,1); end
  if nargout > 3, u = h^2 * (T\f); end

end
