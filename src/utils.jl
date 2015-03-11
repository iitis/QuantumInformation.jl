function funcmh(H::Matrix, f::Function)
  w,v = eig(Hermitian(H))
  fw=diagm(Float64[f(x) for x in w])
  return v*fw*v'
end