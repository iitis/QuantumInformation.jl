random_ginibre_matrix(m::Int,n::Int) = G=randn(n,m)+im*randn(n,m)

function random_dynamical_matrix(n::Int)
  X = random_ginibre_matrix(n^2, n^2)
  Y = ptrace(X*X', [n, n], [1])
  sY = funcmh(Y, x -> 1 / sqrt(x))
  return kron(eye(n,n),sY)*X*X'*kron(eye(n,n),sY)
end