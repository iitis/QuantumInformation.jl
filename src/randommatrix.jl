random_ginibre_matrix!{T<:Float64}(A::Matrix{T}) = randn!(A)

function random_ginibre_matrix!{T<:Complex128}(A::Matrix{T})
    for i=1:length(A)
        A[i] = randn() + 1im * randn()
    end
end

function random_ginibre_matrix{T<:Union(Float64, Complex128)}(M::Type{T}, m::Int64, n::Int64)
    A = zeros(M, m, n)
    random_ginibre_matrix!(A)
    A
end

random_ginibre_matrix(m::Int64, n::Int64) = random_ginibre_matrix(Complex128, m, n)
random_ginibre_matrix(m::Int64) = random_ginibre_matrix(m, m)
random_ginibre_matrix{T<:Union(Float64, Complex128)}(M::Type{T}, m::Int64) = random_ginibre_matrix(M, m, m)

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

function random_dynamical_matrix!{T<:Union(Float64, Complex128)}(J::Matrix{T})
    random_ginibre_matrix!(J)
    n = int(sqrt(size(J, 1)))
    X = A_mul_Bc(J, J)
    Y = ptrace(X, [n, n], [1])
    sY = funcmh!(x -> 1 / sqrt(x), Y)
    A_mul_B!(J, eye(n,n) ⊗ sY * X, eye(n,n) ⊗ sY)
end

function random_dynamical_matrix{T<:Union(Float64, Complex128)}(M::Type{T}, n::Int64)
    J = zeros(M, n*n, n*n)
    random_dynamical_matrix!(J)
    J
end