H = triu(rand(10),-1);
b = rand(10,1);

sigma = rand(1);
xref = (H-sigma*eye(10))\b;
x    = hw7hess(H, sigma, b);

fprintf('Relerr: %e\n', norm(xref-x)/norm(xref));
