dims = 128:128:1024;
dims = sort([dims, dims-1, dims+1]);
ndims = length(dims);
rates = zeros(4,ndims);

for k = 1:length(dims)

  n = dims(k);
  A = rand(n);
  x = rand(n,1);

  ntrials = floor(1e7/n^2);
  gflops = 2*n^2*ntrials/1e9;
  y = A*x;
  y1 = matvec1(A,x);
  y2 = matvec2_row(A,x);
  y3 = matvec2_col(A,x);
  fprintf('%d: %e %e %e\n', n, norm(y-y1), norm(y-y2), norm(y-y3));

  tic;
  for trial = 1:ntrials, y = A*x; end
  rates(1,k) = gflops/toc;

  tic;
  for trial = 1:ntrials, y = matvec1(A,x); end
  rates(2,k) = gflops/toc;

  tic;
  for trial = 1:ntrials, y = matvec2_row(A,x); end
  rates(3,k) = gflops/toc;

  tic;
  for trial = 1:ntrials, y = matvec2_col(A,x); end
  rates(4,k) = gflops/toc;

end

fp = fopen('../data/matvec_time.dat', 'w');
fprintf(fp, 'n matvec0 matvec1 matvec2_row matvec2_col\n');
fprintf(fp, '%d %g %g %g %g\n', [dims; rates]);
fclose(fp);
