random_ginibre_matrix!(A::AbstractMatrix{<:Real}) = randn!(A)

"""
$(SIGNATURES)
- `A`: matrix.

Generates random [Ginibre matrix](https://arxiv.org/pdf/0805.2986.pdf) on the place of matrix `A`.
"""
function random_ginibre_matrix!(A::AbstractMatrix{<:Complex})
    n = prod(size(A))
    A[:] = randn(n) + 1im*randn(n)
end

"""
$(SIGNATURES)
- `m`: dimension.
- `n`: dimension.

Generates a random [Ginibre matrix](https://arxiv.org/pdf/0805.2986.pdf) of dimension `m`x`n`.
"""
function random_ginibre_matrix(::Type{T}, m::Int, n::Int) where T<:Union{Real, Complex}
    A = zeros(T, m, n)
    random_ginibre_matrix!(A)
    A
end

random_ginibre_matrix(m::Int, n::Int) = random_ginibre_matrix(ComplexF64, m, n)
random_ginibre_matrix(m::Int) = random_ginibre_matrix(m, m)
random_ginibre_matrix(::Type{T}, m::Int)  where T<:Union{Real, Complex} = random_ginibre_matrix(T, m, m)

"""
$(SIGNATURES)
- `n`: dimension.

Generates a random `n`x`n` [unitary matrix](https://arxiv.org/pdf/math-ph/0609050.pdf) distributed with Haar measure.
"""
function random_orthogonal(n::Int)
     z = randn(n, n)
     q,r = qr(z)
     d = diag(r)
     ph = d./abs.(d)
     return q.*repmat(ph, 1, size(ph, 1))'
end

"""
$(SIGNATURES)
- `n`: dimension.

Generates a `n`x`n` [Gaussian orthogonal ensemble](https://web.eecs.umich.edu/~rajnrao/Acta05rmt.pdf).
"""
function random_GOE(n::Int)
    H = random_ginibre_matrix(Float64, n)
    return (H+H')/2
end

"""
$(SIGNATURES)
- `n`: dimension.

Generates a `n`x`n` [Gaussian unitary ensemble](https://web.eecs.umich.edu/~rajnrao/Acta05rmt.pdf).
"""
function random_GUE(n::Int)
    H = random_ginibre_matrix(n)
    return (H+H')/2
end

"""
$(SIGNATURES)
- `J`: matrix.
- `A`: matrix.
- `d1`: dimension

Generates random dynamical matrix.
"""
function random_dynamical_matrix!!(J::AbstractMatrix{T}, A::AbstractMatrix{T}, d1::Int=isqrt(size(J, 1))) where T<:Union{Real, Complex}
    random_ginibre_matrix!(A)
    d2 = div(size(J, 1), d1)
    J[:] = A*A'
    Y = ptrace(J, [d2, d1], [1])
    sY = funcmh!(x -> 1 / sqrt(x), Y)
    onesY = eye(d2, d2) ⊗ sY
    J[:] = onesY * J * onesY'
    #can we do this better?
    for i=1:size(J, 1)
        J[i, i] = real(J[i, i])
    end
    J[:] = Hermitian(J)
end

function random_dynamical_matrix!(J::AbstractMatrix{T}, d1::Int=isqrt(size(J,1)), z::Int=prod(size(J))) where T<:Union{Real, Complex}
    A = random_ginibre_matrix(T, size(J, 1), z)
    random_dynamical_matrix!!(J, A)
end

function random_dynamical_matrix(::Type{T}, d1::Int, d2::Int=d1, z::Int=d1*d2) where T<:Union{Real, Complex}
    J = zeros(T, d1*d2, d1*d2)
    random_dynamical_matrix!(J, d1, z)
    J
end

random_dynamical_matrix(d1::Int, d2::Int=d1, z::Int=d1*d2) = random_dynamical_matrix(ComplexF64, d1, d2, z)

"""
$(SIGNATURES)
- `n`: dimension.
- `m`: dimension.

Generates random isometry matrix of dimension `n`x`m`.
"""
function random_isometry(n::Int, m::Int)
    if m > n
        throw(ArgumentError("Isometry must reduce dimensions"))
    end
    z = random_ginibre_matrix(n, m)/sqrt(2.0)
    q,r = qr(z)
    d = diag(r)
    ph = d./abs.(d)
    return q.*repmat(ph, 1, size(q, 1))'
end

"""
$(SIGNATURES)
- `n`: dimension.

Generates random unitary matrix of dimension `n`x`n`.
"""
random_unitary(n::Int) = random_isometry(n, n)
