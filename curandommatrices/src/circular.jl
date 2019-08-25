for T in (ComplexF32, ComplexF64, ComplexF16)
    @eval begin
        @inline CUDAnative.abs(x::$T) = CUDAnative.hypot(x.re, x.im)
    end
end


function _qr_fix!(z::CuMatrix)
    q, r = CuArrays.qr!(z)
    ph = CuArrays.diag(r)
    ph = ph ./ CUDAnative.abs.(ph)
    idim = size(r, 1)
    q = CuMatrix(q)[:, 1:idim]
    q = CuArrays.transpose(ph) .* q
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

function curand(c::CSE)
    z = curand(c.g)
    u = _qr_fix!(z)
    ur = cat([CuMatrix{Float32}([0 -1; 1 0]) for _=1:c.d÷2]..., dims=[1,2])
    ur*u*ur'*transpose(u)
end

for T in (CUE, CircularRealEnsemble, CircularQuaternionEnsemble, HaarIsometry)
    @eval begin
        function curand(c::$T)
            z = curand(c.g)
            _qr_fix!(z)
        end
    end
end