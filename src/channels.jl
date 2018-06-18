
#### Relationship among representations of channels

"""
$(SIGNATURES)
- `kraus`: list of Kraus operators.
- `atol`: tolerance of approximation.

Checks if set of Kraus operators fulfill completness relation.
"""
function kraus_is_CPTP(kraus_list::Vector{<:AbstractMatrix{<:Number}}, atol=1e-08)
    complentess_relation = sum(k'*k for k in kraus_list)
    isapprox(complentess_relation, eye(complentess_relation), atol=atol)
end

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into super-operator matrix.
"""
function kraus_to_superoperator(kraus_list::Vector{<:AbstractMatrix{<:Number}})
    # TODO: chceck if all Kraus operators are the same shape
    sum(k⊗k' for k in kraus_list)
end

"""
$(SIGNATURES)
- `channel`: quantum channel map.
- `dim`: square root of the super-operator matrix dimension.

Transforms quntum channel into super-operator matrix.
"""
function channel_to_superoperator(channel::Function, dim::Int)
    dim > 0 ? () : error("Channel dimension has to be nonnegative")

    M = zeros(ComplexF64, dim*dim, dim*dim)
    for (i, e) in enumerate(base_matrices(dim))
        M[:, i] = res(channel(e))
    end
    M
end

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into Stinespring representation of quantum channel.
"""
function kraus_to_stinespring(kraus_list::Vector{<:AbstractMatrix{<:Number}})
    dim = size(kraus_list,1)
    sum(k⊗ket(i, dim) for (k, i) in enumerate(kraus_list))
end

"""
$(SIGNATURES)
- `kraus_list`: list of Kraus operators.

Transforms list of Kraus operators into dynamical matrix.
"""
function kraus_to_dynamical_matrix(kraus_list::Vector{<:AbstractMatrix{<:Number}})
    sum(res(k) * res(k)' for k in kraus_list)
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into list of Kraus operators.
"""
function superoperator_to_kraus(m::AbstractMatrix{<:Number}, cols=0)

    F = eigfact(Hermitian(reshuffle(m)))
    # TODO: vcat ?
    [sqrt(val)*unres(F.vectors[:,i], cols) for (i,val) in enumerate(F.values)]
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into dynamical matrix.
"""
superoperator_to_dynamical_matrix(m::AbstractMatrix{<:Number}) = reshuffle(m)

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into Stinespring representation of quantum channel.
"""
function superoperator_to_stinespring(m::AbstractMatrix{<:Number})
    kraus_to_stinespring(superoperator_to_kraus(M))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into list of Kraus operators.
"""
function dynamical_matrix_to_kraus(R::AbstractMatrix{<:Number})
    
    F = eigfact(Hermitian(R))

    kraus = Any[] #TODO: make proper type
    for (val, vec) in zip(F.values, F.vectors)
        if val>=0.0
            push!(kraus, sqrt(val) * unres(vec))
        else
            push!(kraus, zero(unres(vec)))
        end
    end
    return kraus
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into Stinespring representation of quantum channel.
"""
function dynamical_matrix_to_stinespring(R::AbstractMatrix{<:Number})
    return kraus_to_stinespring(dynamical_matrix_to_kraus(R))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into super-operator matrix.
"""
function dynamical_matrix_to_superoperator(R::AbstractMatrix{<:Number})
    return reshuffle(R)
end

#= TODO
function general_map_to_kraus(M)
    u, s, v = svd(M)
    s = where(s > 1e-8, s, 0)
    left  = map(x -> sqrt(x[1]) * transpose(unres(x[2])), zip(s, transpose(u))
    right = map(x -> sqrt(x[1]) * unres(conj(x[2])), zip(s, v))
    return [left, right]
end

function general_kraus_to_stinespring(l)
    left, right = l
    num_kraus = len(left)
    basis = [ket(i, num_kraus) for i in range(num_kraus)]
    a0 = sum([kron(left[i], b) for i, b in enumerate(basis)])
    a1 = sum([kron(right[i], b) for i, b in enumerate(basis)])
    return a0, a1
end

function general_map_to_stinespring(R)
    return general_kraus_to_stinespring(general_map_to_kraus(R))
end
=#


#### Application of channels
"""
$(SIGNATURES)
- `R`: dynamical matrix.
- `rho`: quantum state.

Application of dynamical matrix into state `rho`.
"""
function apply_channel_dynamical_matrix(R::AbstractMatrix{<:Number}, rho::AbstractMatrix{<:Number})
    unres(reshuffle(R) * res(rho))
end

"""
$(SIGNATURES)
- `R`: list of Kraus operators.
- `rho`: quantum state.

Application of list of Kraus operators into state `rho`.
"""
function apply_channel_kraus(kraus_list::Vector{<:AbstractMatrix{<:Number}}, rho::AbstractMatrix{<:Number})
    return sum(k * rho * k' for k in kraus_list)
end

"""
$(SIGNATURES)
- `M`: super-operator matrix.
- `rho`: quantum state.

Application of super-operator matrix into state `rho`.
"""
function apply_channel_superoperator(M::AbstractMatrix{<:Number}, rho::AbstractMatrix{<:Number})
    return unres(M * res(rho))
end

"""
$(SIGNATURES)
- `A`: Stinespring representation of quantum channel.
- `rho`: quantum state.
- `dims`: dimensions of registers of `rho`.

Application of Stinespring representation of quantum channel into state `rho`.
"""
function apply_channel_stinespring(A::AbstractMatrix{<:Number}, rho::AbstractMatrix{<:Number}, dims)
    return ptrace(A * rho * A', dims, [1,])
end

function product_superoperator(M1, M2)
    # TODO:
    error("not imeplented")
end
