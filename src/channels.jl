function kraus_check_size(kraus_list::Vector{T}) where T<:AbstractMatrix{<:Number}
    sizes = [size(k) for k in kraus_list]
    for s in sizes[2:end]
        if s!=sizes[1]
            throw(ArgumentError("Kraus operators list contains matrices of different dimmension"))
        end
    end
end
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
    kraus_check_size(kraus_list)
    sum(k⊗(conj.(k)) for k in kraus_list)
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
    dim = size(kraus_list, 1)
    sum(k⊗ket(i-1, dim) for (i, k) in enumerate(kraus_list))
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
function superoperator_to_kraus(M::AbstractMatrix{<:Number})
    F = eigfact(Hermitian(reshuffle(M)))
    _, c = size(M)
    sc = isqrt(c)

    [sqrt(val)*unres(F.vectors[:,i], sc) for (i,val) in enumerate(F.values)]
end

function superoperator_to_kraus(M::AbstractSparseMatrix{<:Number})
    warn("converting to full matrix")
    superoperator_to_kraus(full(M))
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into dynamical matrix.
"""
function superoperator_to_dynamical_matrix(M::AbstractMatrix{<:Number})
    r, c = size(M)
    sr, sc = isqrt(r), isqrt(c)
    reshuffle(M, [sr sr; sc sc])
end

"""
$(SIGNATURES)
- `m`: super-operator matrix.

Transforms super-operator matrix into Stinespring representation of quantum channel.
"""
function superoperator_to_stinespring(M::AbstractMatrix{<:Number})
    kraus_to_stinespring(superoperator_to_kraus(M))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into list of Kraus operators.
"""
function dynamical_matrix_to_kraus(R::Tm, rows, cols) where Tm <:AbstractMatrix{<:Number}
    F = eigfact(Hermitian(R))
    kraus = Tm[]
    for i in 1:length(F.values)
        if F.values[i] >= 0.0
            push!(kraus, sqrt(F.values[i]) * unres(F.vectors[:,i], rows, cols))
        else
            push!(kraus, zero(unres(F.vectors[:,i], rows, cols)))
        end
    end
    return kraus
end

function dynamical_matrix_to_kraus(R::AbstractSparseMatrix{<:Number}, rows, cols)
    warn("converting to full matrix")
    dynamical_matrix_to_kraus(full(R), rows, cols)
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into Stinespring representation of quantum channel.
"""
function dynamical_matrix_to_stinespring(R::AbstractMatrix{<:Number}, rows, cols)
    return kraus_to_stinespring(dynamical_matrix_to_kraus(R, rows, cols))
end

"""
$(SIGNATURES)
- `R`: dynamical matrix.

Transforms dynamical matrix into super-operator matrix.
"""
function dynamical_matrix_to_superoperator(R::AbstractMatrix{<:Number}, rows, cols)
    reshuffle(R, [rows rows; cols cols])
end

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
    return ptrace(A * rho * A', dims, 2)
end

function channel_to_superoperator(channel::Function, dim::Int)
    dim > 0 ? () : error("Channel dimension has to be nonnegative")
    M = zeros(ComplexF64, dim*dim, dim*dim)
    for (i, e) in enumerate(base_matrices(dim))
        M[:, i] = res(channel(e))
    end
    M
end
