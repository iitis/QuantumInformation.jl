export norm_diamond, diamond_distance

const MOI = Convex.MOI

"""

- `Φ`: DynamicalMatrix

Return [diamond norm](https://arxiv.org/pdf/1207.5726.pdf) of dynamical matrix `Φ`.
"""
norm_diamond(Φ::AbstractQuantumOperation, method::Symbol=:primal, eps=1e-7) = norm_diamond(Φ, Val(method), eps)

function norm_diamond(Φ::DynamicalMatrix, ::Val{:primal}, eps)
    J = Φ.matrix
    # TODO: compare d1, d2 with idim, odim
    d1 = Φ.idim
    d2 = Φ.odim
    X = ComplexVariable(d1*d2, d1*d2)
    t = 0.5*inner_product(X, J) + 0.5*inner_product(X', J')

    ρ₀ = ComplexVariable(d1, d1)
    ρ₁ = ComplexVariable(d1, d1)

    constraints = [
        isposdef(ρ₀),
        isposdef(ρ₁),
        tr(ρ₀) == 1,
        tr(ρ₁) == 1,
        isposdef(vcat(hcat(𝕀(d2) ⊗ ρ₀, X), hcat(X', 𝕀(d2) ⊗ ρ₁)))
    ]

    problem = maximize(t, constraints)
    solve!(
        problem,
        MOI.OptimizerWithAttributes(SCS.Optimizer, "eps_abs" => eps);
        silent = true
    )
    problem.optval
end


function norm_diamond(Φ::DynamicalMatrix, ::Val{:dual}, eps)
    J = Φ.matrix
    # TODO: compare d1, d2 with idim, odim
    d1 = Φ.idim
    d2 = Φ.odim
    Y₀ = ComplexVariable(d1*d2, d1*d2)
    Y₁ = ComplexVariable(d1*d2, d1*d2)

	t = 0.5*sigmamax(partialtrace(Y₀, 1, [d2,d1])) +
		0.5*sigmamax(partialtrace(Y₁, 1, [d2,d1]))
    Z = vcat(hcat(Y₀, -J), hcat(-J', Y₁))
    constraints = [
        isposdef(Y₀),
        isposdef(Y₁),
        isposdef(Z+Z')
    ]

    problem = minimize(t, constraints)
    solve!(
        problem,
        MOI.OptimizerWithAttributes(SCS.Optimizer, "eps_abs" => eps);
        silent = true
    )
    problem.optval
end



# function norm_diamond(Φ::AbstractQuantumOperation{T}, args...) where T<:AbstractMatrix{<:Number}
#     norm_diamond(DynamicalMatrix{T}(ϕ), args...)
# end

"""

- `Φ1`: DynamicalMatrix
- `Φ2`: DynamicalMatrix

Return [diamond distance](https://arxiv.org/pdf/1004.4110.pdf) between dynamical matrices `Φ1` and `Φ2`.
"""
function diamond_distance(Φ1::DynamicalMatrix{T}, Φ2::DynamicalMatrix{T}, args...) where T<:AbstractMatrix{<:Number}
    J1 = Φ1.matrix
    J2 = Φ2.matrix
    # TODO: Test dimnesions
    Φ = DynamicalMatrix{T}(J1-J2, Φ1.idim, Φ1.odim)
    norm_diamond(Φ, args...)
end

function diamond_distance(Φ1::AbstractQuantumOperation{T}, Φ2::AbstractQuantumOperation{T}, args...) where T<:AbstractMatrix{<:Number}
    diamond_distance(convert(DynamicalMatrix{T},Φ1), convert(DynamicalMatrix{T},Φ2), args...)
end
