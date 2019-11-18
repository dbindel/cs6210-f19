# Julia code for iterative solves / subspace approx for the model tridiag

using LinearAlgebra
using PyPlot
using FFTW


# Computes a discrete sine transform of the input using FFTW.
#
function dst(f)
  return FFTW.r2r(f, FFTW.RODFT00)
end


# Plots the error in uh as an approximation to u (both on the x grid)
# Also plots the discrete sine transform of the error in order to
# understand it frequency-by-frequency.
#
function plot_err(x, u, uh)
  plt.close()
  subplot(2,1,1)
  plot(x, u, x, uh)
  subplot(2,1,2)
  semilogy(abs.(dst(u[2:end-1]-uh[2:end-1])))
  plt.show(block=false)
end


# Form a basis of cs Chebyshev polynomial evaluations on a uniforrm
# grid of n points.
#
function cheb_basis(n, cs)
  V = zeros(n, cs)
  xx = [1-2*j/(n-1) for j = 0:n-1]
  V[:,1] .= 1
  V[:,2] = xx
  for k = 3:cs
    V[:,k] = 2*(xx.*V[:,k-1]) - V[:,k-2]
  end
  return V
end


# Apply block Gauss-Seidel to update a guess uh for the linear
# system T*u = hf.  bs is the size of the Gauss-Seidel blocks.
#
function solve_bgs(T, uh, hf, bs)
  n = length(uh)-2;
  for j = 2:bs:n+1
    # TODO: Update variables from j:min(j+bs-1.n+1)
  end
  return uh
end


# Use a Bubnov-Galerkin approximation to correct a guess uh
# to the solution of the system T*u = hf, i.e.
# uh += du where T*du = hf-T*u is approximated by Bubnov-Galerkin
# with the basis V.
#
function solve_coarse(T, uh, hf, V)
  r = hf-T*uh;
  # TODO: Apply the Bubnov-Galerkin approximate correction
  return uh
end


# Run tests of the error for a block Gauss-Seidel iteration,
# a Bubnov-Galerkin approximation with a low-dimensional polynomial space,
# and a coarse/fine iteration involving both.
#
function tester()

  # Set up (be lazy and use dense solve)
  n = 100
  h = 1/(n+1)
  x = [j*h for j in 0:n+1]
  f = exp.(5*x)
  hf = h^2*f
  e = ones(n+1)
  T = SymTridiagonal(2*ones(n+2), -ones(n+1))

  # Get a reference solution
  uref = zeros(n+2,1);
  uref[2:n+1] = T[2:n+1,2:n+1]\hf[2:n+1];

  bs = 10;  # Gauss-Seidel block size
  cs = 10;  # Coarse grid space dimension

  # Set up coarse grid problem
  V = cheb_basis(n, cs)
  Tproj = V'*T[2:n+1,2:n+1]*V

  # Plot error for block Gauss-Seidel
  println("-- Block G-S")
  uh = zeros(n+2)
  for sweep = 1:101
    uh = solve_bgs(T, uh, hf, bs)
    if mod(sweep,10) == 1
      plot_err(x, uref, uh);
      println("$(sweep): $(norm(uref-uh))")
      readline()
    end
  end

  # Plot error for a coarse-grid solve
  println("-- Coarse grid");
  uh = zeros(n+2);
  uh = solve_coarse(T, uh, hf, V);
  plot_err(x, uref, uh);

  # TODO: Compute the optimal approximation in the space (uopt)
  #       and the quasi-optimality constant Cproj
  uopt = uh
  Cproj = 1
  
  err = norm(uref-uh)
  err_opt = norm(uref-uopt)
  err_bound = Cproj * err_opt
  println("$(err) $(err_opt) $(err_bound)")
  readline()

  # Plot error for a combined iteration
  println("-- Coarse-fine iteration")
  uh = zeros(n+2)
  for sweep = 1:101
    uh = solve_coarse(T, uh, hf, V)
    uh = solve_bgs(T, uh, hf, bs)
    if mod(sweep,10) == 1
      plot_err(x, uref, uh)
      println("$(sweep): $(norm(uref-uh))")
      readline()
    end
  end
end

tester()