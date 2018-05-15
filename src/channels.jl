
#### Relationship among representations of channels

"""
Checks if set of Kraus operators fulfill completness relation.
"""
function kraus_is_complete(kraus, atol=1e-08)
    complentess_relation = sum(k'*k for k in kraus)
    isapprox(complentess_relation, eye(size(complentess_relation,0)), atol=atol)
end

"""
Transforms list of Kraus operators into super-operator matrix.
"""


function kraus_to_superoperator(kraus_list::Vector{T}) where {T<:AbstractMatrix{T1}} where {T1<:Number}
    # TODO: chceck if all Kraus operators are the same shape
    sum(k⊗k' for k in kraus_list)
end

function channel_to_superoperator(channel::Function, dim::Int)
    dim > 0 ? () : error("Channel dimension has to be nonnegative")

    M = zeros(ComplexF64, dim*dim, dim*dim)
    for (i, e) in enumerate(base_matrices(dim))
        M[:, i] = res(channel(e))
    end
    M
end

function kraus_to_stinespring(kraus::Vector{T}) where {T<:AbstractMatrix{T1}} where {T1<:Number}
    dim = size(kraus,1)
    sum(k⊗ket(i, dim) for (k, i) in enumerate(kraus))
end

function kraus_to_dynamical_matrix(kraus::Vector{T}) where {T<:AbstractMatrix{T1}} where {T1<:Number}
    sum(res(k) * res(k)' for k in kraus)
end

function superoperator_to_kraus(m::T, cols=0) where {T<:AbstractMatrix{T1}} where {T1<:Number}
    F = eigfact(Hermitian(reshuffle(m)))
    # TODO: vcat ?
    [sqrt(val)*unres(F.vectors[:,i], cols) for (i,val) in enumerate(F.values)]
end

superoperator_to_dynamical_matrix(m::T) where {T<:AbstractMatrix{T1}} where {T1<:Number} = reshuffle(m)

function superoperator_to_stinespring(m::T) where {T<:AbstractMatrix{T1}} where {T1<:Number}
    kraus_to_stinespring(superoperator_to_kraus(M))
end

function dynamical_matrix_to_kraus(R):
    vals, vec = scipy.linalg.eigh(R)
    vec_list = [vec[:, i] for i in range(vec.shape[1])]
    filtered = zip(vals, vec_list)  #filter(lambda z: z[0] > 0, zip(vals, vec_list))
    kraus = map(lambda x: np.matrix(np.sqrt(x[0]) * unres(x[1])) if x[0] > 0 else np.matrix(0 * unres(x[1])), filtered)
    return kraus
end

function dynamical_matrix_to_stinespring(R)
    return kraus_to_stinespring(dynamical_matrix_to_kraus(R))
end

function dynamical_matrix_to_superoperator(R)
    return reshuffle(R)
end

function general_map_to_kraus(M)
    u, s, v = svd(M)
    s = np.where(s > 1e-8, s, 0)
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

#### Application of channels

function apply_channel_dynamical_matrix(R, rho)
    unres(reshuffle(R) * res(rho))
end

function apply_channel_kraus(kraus, rho)
    return sum(k * rho * k' for k in kraus)
end

function apply_channel_superoperator(M, rho)
    return unres(M * res(rho))
end

function apply_channel_stinespring(A, rho, dims)
    return ptrace(A * rho * A', dims, [1,])
end

function product_superoperator(M1, M2)
    # TODO:
    error("not imeplented")
end
