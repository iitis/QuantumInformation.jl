function number2mixedradix(n::Int, radices::Vector{Int})
    n >= prod(radices) ? throw(ArgumentError("number to big to transform")) : ()

    digits = Array{Int}(length(radices))
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
    t = trace(ρ)
    for i=1:length(ρ)
        ρ[i] = ρ[i]/t
    end
end

function funcmh!(f::Function, H::Hermitian{T}, R::Matrix{T})  where T<:Union{Real, Complex}
    F = eigfact!(H)
    for i=1:length(F.values)
        F.values[i] = f(F.values[i])
    end
    times_diag = zero(F.vectors)
    for i=1:size(F.vectors, 2)
        times_diag[:, i] = F.vectors[:, i] * F.values[i]
    end
    R[:] = times_diag * F.vectors' # A_mul_Bc!(R, times_diag, F.vectors) deprecated
end

function funcmh!(f::Function, H::Hermitian{T}) where T<:Union{Real, Complex}
    R = zeros(T, size(H))
    funcmh!(f, H, R)
    R
end

function funcmh(f::Function, H::Hermitian{T}) where T<:Union{Real, Complex}
    R = zeros(T, size(H))
    funcmh!(f, copy(H), R)
    R
end

function funcmh!(f::Function, H::Matrix{T}, R::Matrix{T}) where T<:Union{Real, Complex}
    ishermitian(H) ? funcmh!(f, Hermitian(H), R) : error("Non-hermitian matrix passed to funcmh")
end

function funcmh!(f::Function, H::Matrix{T}) where T<:Union{Real, Complex}
    ishermitian(H) ? funcmh!(f, Hermitian(H)) : error("Non-hermitian matrix passed to funcmh")
end

function funcmh(f::Function, H::Matrix{T}) where T<:Union{Real, Complex}
    ishermitian(H) ? funcmh(f, Hermitian(H)) : error("Non-hermitian matrix passed to funcmh")
end

function isidentity(ρ::AbstractMatrix{<:Number}, atol=1e-08)
    rows, cols = size(ρ)
    if rows!=cols
        return false
    end

    isapprox(ρ, eye(ρ), atol=atol)
end

function ispositive(ρ::AbstractMatrix{<:Number}, atol=1e-08)
    rows, cols = size(ρ)
    if rows!=cols
        return false
    end
    if !ishermitian(ρ)
        return false
    end
    fact = eigfact(Hermitian(ρ))
    all(fact.values .> -atol)
end

random_ball(dim::Int) = rand()^(1/dim) * random_ket(Float64, dim)

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
