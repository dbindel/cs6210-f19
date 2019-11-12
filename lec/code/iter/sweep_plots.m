% Set up a model problem
n = 100;
[h,T,f,u] = model2d(n);

% Verify that we get the same behavior with Fourier solver
F = convert_v2m(f);
X = solve_dst(h^2*F(2:end-1,2:end-1));
uu = reshape(X,n*n,1);
fprintf('Test vs FFT solver: %e\n', norm(u-uu)/norm(u));

nsweep = 500;

err_jacobi = zeros(nsweep,1);
U = zeros(n+2,n+2);
for s = 1:nsweep
  U = sweep_jacobi(U,F);
  err_jacobi(s) = norm(U(2:end-1,2:end-1)-X, 'fro');
end

err_gs = zeros(nsweep,1);
U = zeros(n+2,n+2);
for s = 1:nsweep
  U = sweep_gs(U,F);
  err_gs(s) = norm(U(2:end-1,2:end-1)-X, 'fro');
end

err_adi = zeros(nsweep,1);
U = zeros(n+2,n+2);
for s = 1:nsweep
  U = sweep_adi(U,F);
  err_adi(s) = norm(U(2:end-1,2:end-1)-X, 'fro');
end

semilogy(1:nsweep, [err_jacobi, err_gs, err_adi]);
