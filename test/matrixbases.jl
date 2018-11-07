@testset "base_matrices" begin
    d = 4
    m = collect(Matrix{ComplexF64}, base_matrices(4))
    for i=1:d, j=1:d
        v = tr(m[i]' * m[j])
        i == j ? @test(v == 1.) : @test(v == 0.)
    end
end
