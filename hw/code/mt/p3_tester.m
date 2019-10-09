n = 100;
alpha = 1 + rand(n,1);
beta = rand(n-1,1);
T = diag(alpha) + diag(beta,1) + diag(beta,-1);

normTinv_ref = norm(inv(T), 1);
normTinv = p3_norminv(alpha, beta);
relerr = abs(normTinv_ref-normTinv)/normTinv_ref;

printf('3: Relerr %e\n', relerr);
