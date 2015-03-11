include("randommatrix.jl")

function random_mixed_state_hs(d::Int)
  A=random_ginibre_matrix(d,d)
  A=A*A'
  A=A/trace(A)
  return A
end

random_jamiolkowski_state(n::Int) = random_dynamical_matrix(n) / n

function random_ket(d::Int)
  c=randn(d,1)+i*randn(d,1)
  c=c/norm(c)
  return c
end

function random_pure_state(d::Int)
    return proj(random_ket(d))
end