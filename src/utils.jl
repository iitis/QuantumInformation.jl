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
