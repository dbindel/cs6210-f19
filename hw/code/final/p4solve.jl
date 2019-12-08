using LinearAlgebra


# Minimize x'*M*x s.t. Cx = b
#
function p4minx(M, C, b)
end


# Return an M-orthonormal basis V for the null space of C, i.e.
#   V'*M*V = I
#   C*V = 0
#
function p4null(M, C)
end


# Minimize x'*A*X s.t. x'*M*x = 1 and C*x = b.
# Return (x,true) if the problem is feasible, otherwise
#        (xbar,false) with xbar = minimal M-norm solution to C*x=b.
#
function p4optimize(A, M, C, b)
end
