% First test: Check agreement between reference computation and O(n) version

n = 500;
u = rand(n,1);
v = rand(n,1);
x = rand(n,1);

[y1, z1, d1, df1, t1, c1] = hw1p3_ref(u, v, x);
[y2, z2, d2, df2, t2, c2] = hw1p3(u, v, x);

fprintf('Rel differences:\n y:  %e\n z:  %e\n d:  %e\n df: %e\n t:  %e\n c:  %e\n', 
        norm(y1-y2)/norm(y1), norm(z1-z2)/norm(z1),
        norm(d1-d2)/norm(d1), norm(df1-df2)/norm(df1),
        abs(t1-t2)/abs(t1), abs(c1-c2)/abs(c1));

% Second test: Check scalability (probably don't uncomment this until fixed!)

n = 500000;
u = rand(n,1);
v = rand(n,1);
x = rand(n,1);

%[y2, z2, d2, df2, t2, c2] = hw1p3(u, v, x);
