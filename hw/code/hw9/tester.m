% MATLAB ttester for iterative solves / subspace approx for the model tridiag

% Set up (be lazy and use dense solve)
n = 100;
h = 1/(n+1);
x = linspace(0,1,n+2)';
f = exp(5*x);
hf = h^2*f;
e = ones(n+1,1);
T = 2*eye(n+2) - diag(e,-1) - diag(e,1);

% Get a reference solution
uref = zeros(n+2,1);
uref(2:n+1) = T(2:n+1,2:n+1)\hf(2:n+1);

bs = 10;  % Gauss-Seidel block size
cs = 10;  % Coarse grid space dimension

% Set up coarse grid problem
V = cheb_basis(n, cs);
Tproj = V'*T(2:n+1,2:n+1)*V;

% Plot error for block Gauss-Seidel
fprintf('-- Block G-S\n');
uh = zeros(n+2,1);
for sweep = 1:101
  [uh] = solve_bgs(T, uh, hf, bs);
  if mod(sweep,10) == 1
    plot_err(x, uref, uh);
    fprintf('%d: %e\n', sweep, norm(uref-uh));
    pause;
  end
end

% Plot error for a coarse-grid solve
fprintf('-- Coarse grid\n');
uh = zeros(n+2,1);
uh = solve_coarse(T, uh, hf, V);
plot_err(x, uref, uh);

% TODO: Compute uopt = the optimal approximation to u from V
%       Also compute the theoretical quasi-optimality constant
%       (call it Cproj)
uopt = uh;
Cproj = 1;

err = norm(uref-uh);
err_opt = norm(uref-uopt);
fprintf('%e %e %e\n', err, err_opt, Cproj*err_opt);
pause;

% Plot error for a combined iteration
fprintf('-- Coarse-fine iteration\n');
uh = zeros(n+2,1);
for sweep = 1:101
  uh = solve_coarse(T, uh, hf, V);
  uh = solve_bgs(T, uh, hf, bs);
  if mod(sweep,10) == 1
    plot_err(x, uref, uh);
    fprintf('%d: %e\n', sweep, norm(uref-uh));
    pause;
  end
end
