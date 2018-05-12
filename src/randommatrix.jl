random_ginibre_matrix!(A::AbstractMatrix{T}) where T<:Real = randn!(A)

function random_ginibre_matrix!(A::AbstractMatrix{T}) where T<:Complex
    for i=1:length(A)
        A[i] = randn() + 1im * randn()
    end
end

function random_ginibre_matrix(::Type{T}, m::Int, n::Int) where T<:Union{Real, Complex}
    A = zeros(T, m, n)
    random_ginibre_matrix!(A)
    A
end

random_ginibre_matrix(m::Int, n::Int) = random_ginibre_matrix(ComplexF64, m, n)
random_ginibre_matrix(m::Int) = random_ginibre_matrix(m, m)
random_ginibre_matrix(::Type{T}, m::Int)  where T<:Union{Real, Complex} = random_ginibre_matrix(T, m, m)

function random_unitary(n::Int)
  z = random_ginibre_matrix(n,n)/sqrt(2.0)
  q,r = qr(z)
  d = diag(r)
  ph = d./abs.(d)
  return q.*repmat(ph,1,size(ph)[1])'
end

function random_orthogonal(n::Int)
     z = randn(n,n)
     q,r = qr(z)
     d = diag(r)
     ph = d./abs.(d)
     return q.*repmat(ph,1,size(ph)[1])'
end

function random_GOE(n::Int)
    H = randn(n,n)
    return (H+H')/2
end

function random_GUE(n::Int)
    H = randn(n,n)+1im*randn(n,n)
    return (H+H')/2
end

function random_dynamical_matrix!(J::AbstractMatrix{T}) where T<:Union{Real, Complex}
    random_ginibre_matrix!(J)
    n = round(Int, sqrt(size(J, 1)), RoundDown)
    X = J*J' # A_mul_Bc(J, J) deprecated
    Y = ptrace(X, [n, n], [1])
    sY = funcmh!(x -> 1 / sqrt(x), Y)
    onesY = eye(n,n) ⊗ sY
    J[:] = onesY * X * onesY' # A_mul_B!(J, eye(n,n) ⊗ sY * X, eye(n,n) ⊗ sY) deprecated
end

function random_dynamical_matrix(::Type{T}, n::Int) where T<:Union{Real, Complex}
    J = zeros(T, n*n, n*n)
    random_dynamical_matrix!(J)
    J
end

random_dynamical_matrix(n::Int) = random_dynamical_matrix(ComplexF64, n)
