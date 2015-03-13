include("randommatrix.jl")
include("utils.jl")

function random_mixed_state_hs(d::Int)
  A=random_ginibre_matrix(d,d)
  A=A*A'
  A=A/trace(A)
  return A
end

random_jamiolkowski_state(n::Int) = random_dynamical_matrix(n) / n

function random_ket(d::Int)
  c=randn(d,1)+1im*randn(d,1)
  c=c/norm(c)
  return c
end

function random_pure_state(d::Int)
  return proj(random_ket(d))
end

function random_mixed_state_average_purity(d::Int)
  σ = zeros(d)
  while true
      ρ = random_mixed_state_hs(d)
      purity = trace(ρ * ρ)
      p = real(sqrt(d^2 - 1) / (sqrt(d^2 + 1) * sqrt(d * purity - 1)))
      if p >= 0 && p <= 1
          σ = p * ρ + (1 - p) / d * eye(d)
          break
      end
  end
  return σ
end

function random_mixed_state_fixed_purity(d::Int, p::Real)
  error("Not implemented")
  l = random_vector_fixed_l1_l2(1., p, d)
  u = random_unitary(d)
  return u' * diagm(l) * u
end