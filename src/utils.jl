function renormalize!{T<:Union(Float64, Complex128)}(ϕ::Vector{T})
    n = norm(ϕ)
    for i=1:length(ϕ)
        ϕ[i] = ϕ[i]/n
    end
end

function renormalize!{T<:Union(Float64, Complex128)}(ρ::Matrix{T})
    t = trace(ρ)
    for i=1:length(ρ)
        ρ[i] = ρ[i]/t
    end
end

function funcmh(H::Matrix, f::Function)
  w,v = eig(Hermitian(H))
  fw=diagm(Float64[f(x) for x in w])
  return v*fw*v'
end

function random_sphere(dim::Int)
  v = randn(dim)
  return v / norm(v)
end

random_ball(dim::Int) = rand()^(1/dim) * random_sphere(dim)

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