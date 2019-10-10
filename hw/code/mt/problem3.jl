using LinearAlgebra

# Compute the norm(inv(T), 1) where T is an spd tridiagonal
# with diagonal alpha and off-diagonal beta
#
function norminv(alpha, beta)
  # TODO: Replace with something efficient (can do in O(n))
  return norm(inv(SymTridiagonal(alpha, beta)), 1)
end


# -- Test harness

n = 100
alpha = 1 .+ rand(Float64, n)
beta = rand(Float64, n-1)
T = SymTridiagonal(alpha, beta)

normTinv_ref = norm(inv(T), 1)
normTinv = norminv(alpha, beta)
relerr = abs(normTinv_ref-normTinv)/normTinv_ref

println("3: Relerr $(relerr)")