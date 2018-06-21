struct DynamicalMat{T}
    matrix::T
end

struct Superoperator{T}
    matrix::T
end

struct KrausOps{T}
    matrices::Vector{T}
end

function KrausOps(kraus_list::Vector{T}) where T<:AbstractMatrix{<:Number}
    sizes = [size(k) for k in kraus_list]

    KrausOps
end

#### Relationship among representations of channels

"""
Checks if set of Kraus operators fulfill completness relation.
"""
function kraus_is_CPTP(kraus_list::Vector{<:AbstractMatrix{<:Number}}, atol=1e-08)
    complentess_relation = sum(k'*k for k in kraus_list)
    isapprox(complentess_relation, eye(complentess_relation), atol=atol)
end

"""
Transforms list of Kraus operators into super-operator matrix.
"""
function kraus_to_superoperator(kraus_list::Vector{<:AbstractMatrix{<:Number}})
    # TODO: chceck if all Kraus operators are the same shape
    sum(k⊗(conj.(k)) for k in kraus_list)
end

function kraus_to_stinespring(kraus_list::Vector{<:AbstractMatrix{<:Number}})
    dim = size(kraus_list, 1)
    sum(k⊗ket(i-1, dim) for (i, k) in enumerate(kraus_list))
end

function kraus_to_dynamical_matrix(kraus_list::Vector{<:AbstractMatrix{<:Number}})
    sum(res(k) * res(k)' for k in kraus_list)
end

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

function superoperator_to_dynamical_matrix(M::AbstractMatrix{<:Number})
    r, c = size(M)
    sr, sc = isqrt(r), isqrt(c)
    reshuffle(M, [sr sr; sc sc])
end

function superoperator_to_stinespring(M::AbstractMatrix{<:Number})
    # TODO: This is no the right way for transoformations
    kraus_to_stinespring(superoperator_to_kraus(M))
end

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

function dynamical_matrix_to_kraus(R::AbstractSparseMatrix{<:Number})
    warn("converting to full matrix")
    dynamical_matrix_to_kraus(full(M))
end

function dynamical_matrix_to_stinespring(R::AbstractMatrix{<:Number})
    # TODO: This is no the right way for transoformations
    return kraus_to_stinespring(dynamical_matrix_to_kraus(R))
end

function dynamical_matrix_to_superoperator(R::AbstractMatrix{<:Number})
    r, c = size(R)
    sr, sc = isqrt(r), isqrt(c)
    reshuffle(R, [sr sc; sr sc]) # ???
end

#### Application of channels

function apply_channel_kraus(kraus_list::Vector{<:AbstractMatrix{<:Number}}, rho::AbstractMatrix{<:Number})
    return sum(k * rho * k' for k in kraus_list)
end

function apply_channel_superoperator(M::AbstractMatrix{<:Number}, rho::AbstractMatrix{<:Number})
    return unres(M * res(rho))
end

function apply_channel_stinespring(A::AbstractMatrix{<:Number}, rho::AbstractMatrix{<:Number}, dims)
    return ptrace(A * rho * A', dims, 1)
end

function apply_channel_dynamical_matrix(R::AbstractMatrix{<:Number}, rho::AbstractMatrix{<:Number})
    unres(reshuffle(R) * res(rho))
end

function channel_to_superoperator(channel::Function, dim::Int)
    dim > 0 ? () : error("Channel dimension has to be nonnegative")
    M = zeros(ComplexF64, dim*dim, dim*dim)
    for (i, e) in enumerate(base_matrices(dim))
        M[:, i] = res(channel(e))
    end
    M
end
