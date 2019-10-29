H0 = triu(rand(10),-1);
u = rand(10,1);
v = rand(10,1);
A = H0+u*v';

[Q, H] = hw7hess(H0, u, v);
err1   = norm(A-Q*H*Q')/norm(A);
err2   = norm(tril(H,-2));
err3   = norm(Q'*Q-eye(10));

fprintf('Check recovery of A: %e\n', err1);
fprintf('Check Hessenberg H:  %e\n', err2);
fprintf('Check orthogonal Q:  %e\n', err3);
