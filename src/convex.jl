using Convex, SCS

function diamond_norm(J::AbstractMatrix{<:Number}, d1::Int, d2::Int=d1)
  X = ComplexVariable(d1*d2, d1*d2)
  t = 0.5*inner_product(X, J) + 0.5*inner_product(X', J')

  ρ₀ = ComplexVariable(d1, d1)
  ρ₁ = ComplexVariable(d1, d1)

  constraints = [ρ₀ in :SDP, ρ₁ in :SDP]
  constraints += trace(ρ₀) == 1
  constraints += trace(ρ₁) == 1
  constraints += [eye(d2) ⊗ ρ₀ X; X' eye(d2) ⊗ ρ₁] in :SDP

  problem = maximize(t, constraints)
  solve!(problem, SCSSolver(verbose=0))
  problem.optval
end

function diamond_distance(J1::AbstractMatrix{<:Number}, J2::AbstractMatrix{<:Number}, d1::Int, d2::Int=d1)
    diamond_norm(J1-J2, d1, d2)
end
