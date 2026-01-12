export number2mixedradix, mixedradix2number, funcmh, funcmh!, renormalize!
# realdiag, realdiag!
using SparseArrays

"""
  - `n`: Number to be converted (integer).
  - `radices`: Vector of mixed radices.

Returns the representation of `n` in the mixed radix system defined by `radices`.
"""
function number2mixedradix(n::Int, radices::Vector{Int})
    n >= prod(radices) ? throw(ArgumentError("number to big to transform")) : ()

    digits = Array{Int}(undef, length(radices))
    for (i, radix) in enumerate(reverse(radices))
        n, digits[end - i + 1] = divrem(n, radix)
    end
    return digits
end

"""
  - `digits`: Vector of coefficients in mixed radix representation.
  - `radices`: Vector of mixed radices.

Returns the integer number corresponding to the mixed radix representation.
"""
function mixedradix2number(digits::Vector{Int}, radices::Vector{Int})
    length(digits)>length(radices) ? throw(ArgumentError("more digits than radices")) : ()

    res = 0
    digitsreversed = reverse(digits)
    for (digit, radix) in zip(digits, radices)
        digit >= radix ? throw(ArgumentError("digit larger or equal to base")) : ()
        res = res * radix + digit
    end
    return res
end

"""
  - `ψ`: Input vector.

Renormalizes the vector `ψ` in-place so that its norm is 1.
"""
function renormalize!(ψ::AbstractVector{<:Number})
    n = norm(ψ)
    for i in 1:length(ψ)
        ψ[i] = ψ[i]/n
    end
end

"""
  - `ρ`: Input matrix.

Renormalizes the matrix `ρ` in-place so that its trace is 1.
"""
function renormalize!(ρ::AbstractMatrix{<:Number})
    t = tr(ρ)
    for i in 1:length(ρ)
        ρ[i] = ρ[i]/t
    end
end

function renormalize!(ψ::AbstractSparseVector{<:Number})
    n = norm(ψ)
    nz = nonzeros(ψ)
    for i in 1:length(nz)
        nz[i] /= n
    end
end

function renormalize!(ρ::AbstractSparseMatrix{<:Number})
    t = tr(ρ)
    nz = nonzeros(ρ)
    for i in 1:length(nz)
        nz[i] /= t
    end
end

#FIXME: here be dragons again
# function realdiag!(a::AbstractMatrix{ComplexF64})
#     r, c = size(a)
#     r == c ? () : throw(ArgumentError("Non-square matrix"))
#     for i=1:r
#         a[i, i] = real(a[i, i])
#     end
# end
#
# function realdiag(a::AbstractMatrix{ComplexF64})
#     b = copy(a)
#     realdiag!(b)
#     b
# end
#
# function realdiag(a::AbstractMatrix{<:Number})
#     a
# end

function funcmh!(
    f::Function,
    h::Hermitian{T},
    r::Matrix{T},
) where {T <: Union{Real, Complex}}
    fact = eigen!(h)
    times_diag = zero(fact.vectors)
    for i in 1:size(fact.vectors, 2)
        times_diag[:, i] = fact.vectors[:, i] * f(fact.values[i])
    end
    return r[:] = times_diag * fact.vectors'
end

function funcmh!(f::Function, h::Hermitian{T}) where {T <: Union{Real, Complex}}
    r = zeros(T, size(h))
    funcmh!(f, h, r)
    return r
end

function funcmh(f::Function, h::Hermitian{T}) where {T <: Union{Real, Complex}}
    r = zeros(T, size(h))
    funcmh!(f, copy(h), r)
    return r
end

function funcmh!(f::Function, h::Matrix{T}, r::Matrix{T}) where {T <: Union{Real, Complex}}
    return ishermitian(h) ? funcmh!(f, Hermitian(h), r) :
           error("Non-hermitian matrix passed to funcmh")
end

function funcmh!(f::Function, h::Matrix{T}) where {T <: Union{Real, Complex}}
    return ishermitian(h) ? funcmh!(f, Hermitian(h)) :
           error("Non-hermitian matrix passed to funcmh")
end

function funcmh(f::Function, h::Matrix{T}) where {T <: Union{Real, Complex}}
    return ishermitian(h) ? funcmh(f, Hermitian(h)) :
           error("Non-hermitian matrix passed to funcmh")
end

"""
  - `ρ`: Input matrix.
  - `atol`: Absolute tolerance.

Checks if the matrix `ρ` is approximately the identity matrix.
"""
function isidentity(ρ::AbstractMatrix{<:Number}, atol=1e-13)
    rows, cols = size(ρ)
    if rows!=cols
        return false
    end

    return isapprox(ρ, I, atol=atol)
end

"""
  - `ρ`: Input matrix.
  - `atol`: Absolute tolerance.

Checks if the matrix `ρ` is positive semi-definite.
"""
function ispositive(ρ::AbstractMatrix{<:Number}, atol=1e-13)
    rows, cols = size(ρ)
    if rows!=cols
        return false
    end

    if issparse(ρ)
        ρ = Array(ρ)
    end

    # if !ishermitian(ρ) # TODO: ishermitian function has no tolerance
    #     return false
    # end
    h = Hermitian(ρ)
    fact = eigen(h)
    return all(fact.values .> -atol)
end

isnumbernotint(T::Type) = ((T <: Real && !(T <: Integer)) || (T <: Complex))

#function random_vector_fixed_l1_l2(l1::Real, l2::Real, d::Int)
#  #from here http://stats.stackexchange.com/questions/61692/generating-vectors-under-constraints-on-1-and-2-norm
#  u, _ = qr(ones(d, d))
#  u = -u
#  z = random_sphere(d - 1)
#  z = [0; z]
#  r = sqrt(l2 - l1^2 / d)
#  v = u * z * r
#  return v + l1 / d * ones(d)
#end
