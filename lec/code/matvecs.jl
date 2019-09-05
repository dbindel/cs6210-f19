# Different matvec organizations

# We need LinearAlgebra for the norm function, and Printf for printf.
using LinearAlgebra
using Printf


"""
    matvec_row(A, x)

Compute a matrix-vector product via dot products (i, j order).
"""
function matvec_row(A, x)
  m, n = size(A)
  y = zeros(eltype(x), m)
  for i = 1:m
    for j = 1:n
      y[i] += A[i,j] * x[j]
    end
  end
  return y
end



"""
    matvec_col(A, x)

Compute a matrix-vector product as a sum of scaled columns (j, i order)
"""
function matvec_col(A, x)
  m, n = size(A)
  y = zeros(eltype(x), m)
  for j = 1:n
    for i = 1:m
      y[i] += A[i,j] * x[j]
    end
  end
  return y
end


"""
    matvec_timing(n)

Time the built-in matvec against matvec_row and matvec_col for a random
problem of size n.  Print a measure of the difference and return a string
with the timing information.
"""
function matvec_timing(n)

  A = rand(n,n)
  x = rand(n)
  
  ntrials = floor(1e7 / n^2)
  gflops = 2 * n^2 * ntrials / 1e9
  
  y = A*x
  yr = matvec_row(A, x)
  yc = matvec_col(A, x)
  @printf("%d: %e %e\n", n, norm(y-yr), norm(y-yc))

  t1 = @elapsed begin
    for trial = 1:ntrials
      y = A*x
    end
  end

  t2 = @elapsed begin
    for trial = 1:ntrials
      y = matvec_row(A,x)
    end
  end

  t3 = @elapsed begin
    for trial = 1:ntrials
      y = matvec_col(A,x)
    end
  end

  return @sprintf("%d %g %g %g", n, gflops/t1, gflops/t2, gflops/t3);
end


# Write a data file with timings for a range of sizes
open("../data/matvec_time_jl.dat", "w") do f
  for d = 128:128:1024
    println(f, matvec_timing(d-1))
    println(f, matvec_timing(d))
    println(f, matvec_timing(d+1))
  end
end