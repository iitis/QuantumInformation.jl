using CuArrays, CuTensorOperations, QuantumInformation

for (fname, f) in ((:bench_device, QuantumInformation.curand), (:bench_cpu, rand))
    @eval begin
        function $fname(d=256)
            c = HaarKet(d^2)
            psi = $f(c)
            @time ptrace(psi, [d, d], 1)
            @time ptrace(psi, [d, d], 2)
        end
    end
end

bench_device()
bench_device()
bench_cpu()
bench_cpu()