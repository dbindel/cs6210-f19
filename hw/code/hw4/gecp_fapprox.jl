using LinearAlgebra


# Approximate a function f on [-1,1]^2 via (at most) m steps of
# Gaussian elimination with complete pivoting.
#
# Arguments:
#     f: Function to be approximated
#     m: Maximum number of terms (default 20)
#     rtol (keyword, 1e-8): Relative tolerance for termination
#     atol (keyword, 1e-8): Absolute tolerance for termination
#
# Returns:
#     X: array of X coordinates of pivot points
#     Y: array of Y coordinates of pivot points
#     fpiv: values of the error function (at the points of max modulus)
#
function gecp_fapprox(f, m=20; rtol=1e-8, atol=1e-8)

    # DSB: In my implementation, I maintained an LU factorization
    #      of the matrix with entries f(X[i], Y[j]), stored in packed
    #      form.  I also kept the G and H matrices, computed one column/row
    #      at a time.  I only searched for X and Y over a grid of points;
    #      something more robust might start with a grid search and then
    #      refine with Newton iteration (for example).
 
    fpiv = zeros(Float64, m)
    X = zeros(Float64, m)
    Y = zeros(Float64, m)
    LU_fXY = zeros(Float64, m, m)

    # Set up sample mesh and evaluate f on grid
    nmesh = 101
    xmesh = range(-1, 1, length=nmesh)

    S = [f(x,y) for x=xmesh, y=xmesh]
    G = zeros(Float64, nmesh, m)
    H = zeros(Float64, m, nmesh)
    
    # GECP main loop
    for k = 1:m

        # TODO: Update X, Y, LU_fXY, G, H, and S; quit early if
        #       tolerances are satisfied.
 
    end
    return X, Y, fpiv

end


# Print diagnostics
#
function print_diagnostics(X, Y, fpiv)
    for k = 1:length(X)
        print("$k: err($(X[k]), $(Y[k])) = $(fpiv[k])\n")
    end
end


# === Test things out on some sample functions ===

print("\n=== Squared exp ===\n")
bump(x, y) = exp(-(x^2+y^2)/2)
X, Y, fpiv = gecp_fapprox(bump, 20)
print_diagnostics(X, Y, fpiv)

print("\n=== Ackley ===\n")
ackley(x, y) = 20-20*exp(-0.2*sqrt((x^2+y^2)/2))-exp((cos(2*pi*x)+cos(2*pi*y))/2)
X, Y, fpiv = gecp_fapprox(ackley, 20)
print_diagnostics(X, Y, fpiv)
