random_ginibre_matrix(m::Int,n::Int) = G=randn(n,m)+im*randn(n,m)

function random_unitary(n::Integer)
  z = random_ginibre_matrix(n,n)/sqrt(2.0)
  q,r = qr(complex128(z))
  d = diag(r)
  ph = d./abs(d)
  return q.*repmat(ph,1,size(ph)[1])'
end

function random_orthogonal(n::Integer)
     z = randn(n,n)
     q,r = qr(z)
     d = diag(r)
     ph = d./abs(d)
     return q.*repmat(ph,1,size(ph)[1])'
end

function random_GOE(n::Integer)
    H=randn(n,n)
    return (H+H')/2
end

function random_GUE(n::Integer)
    H=randn(n,n)+1im*randn(n,n)
    return (H+H')/2
end

function random_dynamical_matrix(n::Int)
  X = random_ginibre_matrix(n^2, n^2)
  Y = ptrace(X*X', [n, n], [1])
  sY = funcmh(Y, x -> 1 / sqrt(x))
  return kron(eye(n,n),sY)*X*X'*kron(eye(n,n),sY)
end