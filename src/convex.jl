using Convex, SCS

"""
$(SIGNATURES)
- `J`: matrix.
- `d1`: ?.
- `d2`: ?.

Return [diamond norm](https://arxiv.org/pdf/1004.4110.pdf) of matrix `J`.
"""

function diamond_norm(Φ::AbstractQuantumOperation{T}) where T<:AbstractMatrix{<:Number}
    J = DynamicalMatrix{T}(Φ).matrix
    # TODO: compare d1, d2 with idim, odim
    d1 = Φ.idim
    d2 = Φ.odim
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

"""
$(SIGNATURES)
- `J1`: matrix.
- `J2`: matrix.
- `d1`: ?.
- `d2`: ?.

Return [diamond distance](https://arxiv.org/pdf/1004.4110.pdf) between matrices `J1` and `J2`.
"""
function diamond_distance(Φ1::AbstractQuantumOperation{T}, Φ2::AbstractQuantumOperation{T}) where T<:AbstractMatrix{<:Number}
    J1 = DynamicalMatrix{T}(Φ1).matrix
    J2 = DynamicalMatrix{T}(Φ2).matrix
    # TODO: Test dimnesions
    Φ = DynamicalMatrix{T}(J1-J2, Φ1.idim, Φ1.odim)
    diamond_norm(Φ)
end
