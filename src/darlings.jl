# Functions in this file shall not be used.  They are untested and probably do
# not work as intended


## Channels

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

function product_superoperator(M1, M2)
    # TODO:
    error("not imeplented")
end
