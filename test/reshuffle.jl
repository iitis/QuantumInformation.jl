@testset "Reshuffle" begin
@testset "Dense matrices" begin
    X = reshape([1:16;], 4, 4)'
    T = [1 2 5 6; 3 4 7 8; 9 10 13 14; 11 12 15 16]
    @test reshuffle(X) == T

    X = reshape(1:36, 6, 6)'
    T=[1	2	7	8	13	14;
    3	4	9	10	15	16;
    5	6	11	12	17	18;
    19	20	25	26	31	32;
    21	22	27	28	33	34;
    23	24	29	30	35	36]

    T2 = [1	2	3	7	8	9	13	14	15;
    4	5	6	10	11	12	16	17	18;
    19	20	21	25	26	27	31	32	33;
    22	23	24	28	29	30	34	35	36]

    @test reshuffle(X, [2 3; 3 2]) == T
    @test reshuffle(X, [2 3; 2 3]) == T2

end

end
