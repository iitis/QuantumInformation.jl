export number2mixedradix, mixedradix2number,
    funcmh, funcmh!, renormalize!
    # realdiag, realdiag!

function number2mixedradix(n::Int, radices::Vector{Int})
    n >= prod(radices) ? throw(ArgumentError("number to big to transform")) : ()

    digits = Array{Int}(undef, length(radices))
    for (i, radix) in enumerate(reverse(radices))
        n, digits[end-i+1] = divrem(n, radix)
    end
    digits
end

function mixedradix2number(digits::Vector{Int}, radices::Vector{Int})
    length(digits)>length(radices) ? throw(ArgumentError("more digits than radices")) : ()

    res = 0
    digitsreversed = reverse(digits)
    for (digit, radix) = zip(digits, radices)
        digit >= radix ? throw(ArgumentError("digit larger or equal to base")) : ()
        res = res * radix + digit
    end
    res
end

function renormalize!(ψ::AbstractVector{<:Number})
    n = norm(ψ)
    for i=1:length(ψ)
        ψ[i] = ψ[i]/n
    end
end

function renormalize!(ρ::AbstractMatrix{<:Number})
    t = tr(ρ)
    for i=1:length(ρ)
        ρ[i] = ρ[i]/t
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

function funcmh!(f::Function, h::Hermitian{T}, r::Matrix{T})  where T<:Union{Real, Complex}
    fact = eigen!(h)
    times_diag = zero(fact.vectors)
    for i=1:size(fact.vectors, 2)
        times_diag[:, i] = fact.vectors[:, i] * f(fact.values[i])
    end
    r[:] = times_diag * fact.vectors'
end

function funcmh!(f::Function, h::Hermitian{T}) where T<:Union{Real, Complex}
    r = zeros(T, size(h))
    funcmh!(f, h, r)
    r
end

function funcmh(f::Function, h::Hermitian{T}) where T<:Union{Real, Complex}
    r = zeros(T, size(h))
    funcmh!(f, copy(h), r)
    r
end

function funcmh!(f::Function, h::Matrix{T}, r::Matrix{T}) where T<:Union{Real, Complex}
    ishermitian(h) ? funcmh!(f, Hermitian(h), r) : error("Non-hermitian matrix passed to funcmh")
end

function funcmh!(f::Function, h::Matrix{T}) where T<:Union{Real, Complex}
    ishermitian(h) ? funcmh!(f, Hermitian(h)) : error("Non-hermitian matrix passed to funcmh")
end

function funcmh(f::Function, h::Matrix{T}) where T<:Union{Real, Complex}
    ishermitian(h) ? funcmh(f, Hermitian(h)) : error("Non-hermitian matrix passed to funcmh")
end

function isidentity(ρ::AbstractMatrix{<:Number}; atol=1e-13)
    rows, cols = size(ρ)
    if rows!=cols
        return false
    end

    isapprox(ρ, I, atol=atol)
end

function ispositive(ρ::AbstractMatrix{<:Number}; atol=1e-13)
    rows, cols = size(ρ)
    if rows!=cols
        return false
    end
    # if !ishermitian(ρ) # TODO: ishermitian function has no tolerance
    #     return false
    # end
    h = Hermitian(ρ)
    fact = eigen(h)
    all(fact.values .> -atol)
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
