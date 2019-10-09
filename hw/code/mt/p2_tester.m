
% -- Finite difference check for dnormx

x = rand(10, 1);
A = rand(10, 10);
M = A'*A;
dx = rand(10, 1);
dM = rand(10, 10);
dM = (dM + dM')/2;

h = 1e-6;
normx_p = (x+h*dx)'*(M+h*dM)*(x+h*dx);
normx_m = (x-h*dx)'*(M-h*dM)*(x-h*dx);
dnormx_fd = (normx_p - normx_m)/(2*h);

dnormx = p2_deriv_mnorm(x, M, dx, dM);
relerr = abs(dnormx_fd-dnormx)/abs(dnormx);
printf('2a: Relerr (vs finite diff): %e\n', relerr);


% -- Finite difference check for solve

n = 100;
k = 5;
Z = rand(n, k);
dZ = rand(n, k);
[Q, R] = qr(Z, 0);
b = rand(n, 1);
db = rand(n, 1);

Zp = Z+h*dZ;
Zm = Z-h*dZ;
bp = b+h*db;
bm = b-h*db;
x =  (eye(n) + Z*Z')\b;
xp = (eye(n) + Zp*Zp')\bp;
xm = (eye(n) + Zm*Zm')\bm;
dx_fd = (xp - xm)/(2*h);

dx = p2_deriv_solve(Q, R, x, b, dZ, db);
relerr = norm(dx_fd-dx)/norm(dx);
printf('2b: Relerr (vs finite diff): %e\n', relerr);


% -- Finite difference check for Cholesky

n = 10;
A = rand(n, n);
M = A'*A;
R = chol(M, 'upper');
dM = rand(n, n);
dM = (dM + dM')/2;

Rp = chol(M+h*dM, 'upper');
Rm = chol(M-h*dM, 'upper');
dR_fd = (Rp-Rm)/(2*h);

dR = p2_deriv_chol(M, R, dM);
relerr = norm(dR_fd-dR)/norm(dR);
printf('2c: Relerr (vs finite diff): %e\n', relerr);
