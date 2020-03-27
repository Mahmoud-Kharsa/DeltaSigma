using MAT

@testset "scaleABCD" begin
    # MATLAB output, generated using output of Julia randn
    test = matread("resources/scaleABCD.mat")

    H = synthesizeNTF(5, 42, 1)
    a, g, b, c = realizeNTF(H)
    ABCD = stuffABCD(a, g, b, c)
    ABCDs, umax, S = scaleABCD(ABCD, 2, 0, 1, NaN, NaN, 10000)
    @test isapprox(test["ABCDs"], ABCDs)
    @test isapprox(test["umax"], umax)
    @test isapprox(test["S"], S)
end
