A = rand(10);
b = rand(10,1);
c = rand(10,1);
sigma = 5;

[sigma, gs] = hw7newton(A, b, c, sigma);
err = min(abs(eig(A)-sigma));
fprintf('Final error: %e\n', err);
fprintf('Values for g:\n');
fprintf('  %e\n', gs);
