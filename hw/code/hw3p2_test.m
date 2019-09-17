n = 10;
d = rand(n,1);
b = rand(n,1);
c = rand(n,1);
f = rand(1);
y = rand(n+1,1);
[detA_ref, x_ref] = hw3p2_ref(d, b, c, f, y);
[detA, x] = hw3p2(d, b, c, f, y);

fprintf('relerr detA: %e\n', abs(detA_ref-detA)/abs(detA));
fprintf('relerr x:    %e\n', norm(x-x_ref)/norm(x));
