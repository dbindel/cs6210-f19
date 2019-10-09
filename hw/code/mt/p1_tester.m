n = 100;
k = 5;
Z = rand(n, k);
A = eye(n) + Z*Z';

x = rand(n, 1);
yref = A*x;
y = p1_fast_mul(Z, x);
relerr = norm(y-yref)/norm(yref);
printf('1a: Relative err %e\n', relerr)

[U, S, V] = svd(Z, 0);
kappa_ref = cond(A, 2);
kappa = p1_fast_cond(U, S, V);
relerr = abs(kappa-kappa_ref)/kappa_ref;
printf('1b: Relative err %e\n', relerr)

[Q, R] = qr(Z, 0);
b = rand(n, 1);
xref = A\b;
x = p1_fast_solve(Q, R, b);
relerr = norm(x-xref)/norm(xref);
printf('1c: Relative err %e\n', relerr);
