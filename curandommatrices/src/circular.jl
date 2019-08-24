# FIXME: this can be accelerated
function cplx_phase!(a)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    @inbounds a[i] = a[i] / sqrt(real(a[i])^2 + imag(a[i])^2)
    return nothing
end

# FIXME: use blocks to support larger than 1024

function _qr_fix!(z::CuMatrix)
    q, r = CuArrays.qr!(z)
    ph = diag(r)
    len = length(ph)
    @cuda threads=len cplx_phase!(ph)
    q = CuMatrix(q)
    idim = size(r, 1)
    for i=1:idim
        q[:, i] .*= ph[i]
    end
    q[:, 1:idim]
end

function _qr_fix(z::CuMatrix)
    a = copy(z)
    _qr_fix!(a)
end

function curand(c::COE)
    z = curand(c.g)
    u = _qr_fix!(z)
    transpose(u)*u
end

function curand(c::CUE)
    z = curand(c.g)
    u = _qr_fix!(z)
    u
end

function curand(c::CSE)
    z = curand(c.g)
    u = _qr_fix!(z)
    ur = cat([CuMatrix{Float32}([0 -1; 1 0]) for _=1:c.d÷2]..., dims=[1,2])
    ur*u*ur'*transpose(u)
end

# struct CircularRealEnsemble <: QIContinuousMatrixDistribution
#     d::Int
#     g::GinibreEnsemble{1}

#     function CircularRealEnsemble(d::Int)
#         g = GinibreEnsemble{1}(d)
#         new(d, g)
#     end
# end

# function rand(rng::AbstractRNG, c::CircularRealEnsemble)
#     z = rand(rng, c.g)
#     _qr_fix!(z)
# end

# struct CircularQuaternionEnsemble <: QIContinuousMatrixDistribution
#     d::Int
#     g::GinibreEnsemble{4}

#     function CircularQuaternionEnsemble(d::Int)
#         g = GinibreEnsemble{4}(d)
#         new(d, g)
#     end
# end

# function rand(rng::AbstractRNG, c::CircularQuaternionEnsemble)
#     z = rand(rng, c.g)
#     _qr_fix!(z)
# end


# struct HaarIsometry <: QIContinuousMatrixDistribution
#     idim::Int
#     odim::Int
#     g::GinibreEnsemble{2}

#     function HaarIsometry(idim::Int, odim::Int)
#         idim <= odim || throw(ArgumentError("idim can't be greater than odim"))
#         g = GinibreEnsemble{2}(odim, idim)
#         new(idim, odim, g)
#     end
# end

# function rand(rng::AbstractRNG, c::HaarIsometry)
#     z = rand(rng, c.g)
#     _qr_fix!(z)
# end