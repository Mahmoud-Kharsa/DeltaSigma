using ControlSystems
using MAT

@testset "simulateMS" begin
    # MATLAB output
    test = matread("resources/simulateMS.mat")

    # input parameters
    M = 16
    N = 16384
    f = 164
    A0 = undbv(-3)
    A1 = undbv(-30)
    ntf = synthesizeNTF(6, 25, 1, 4)
    mtf1 = zpk([1], [0], 1, 1)
    mtf2 = zpk([1, 1], [0.3, 0.3], 1, 1)
    mtf3 = synthesizeNTF(2, 25, 1, 2)
    mtf4 = synthesizeNTF(4, 22.5, 1, 1.3)
    A   = [  A0,   A0,  A1,    A0,   A0]
    mtf = [mtf1, mtf2, mtf1, mtf3, mtf4]

    for i = 1:length(mtf)
        u = M * A[i] * sin.(2*pi*f/N * (0:N-1))
        v, = simulateDSM(u, ntf, M+1)
        sv, sx, sigma_se, max_sx, max_sy = simulateMS(v, M, mtf[i])
        @test test["sv_$i"] == sv
        @test isapprox(test["sx_$i"], sx)
        @test isapprox(test["sigma_se_$i"], sigma_se)
        @test isapprox(test["max_sx_$i"], max_sx)
        @test isapprox(test["max_sy_$i"], max_sy)
    end

    # test with dw != [1, 1, ... , 1]
    dw = [8; 8; 4; 4; 2; 2; 1; 1; zeros(Int, 8)]
    u = M * A0 * sin.(2*pi*f/N * (0:N-1))
    v, = simulateDSM(u, ntf, M+1)
    sv, sx, sigma_se, max_sx, max_sy = simulateMS(v, M, mtf4, 0 , dw)
    @test test["sv_dw"] == sv
    @test isapprox(test["sx_dw"], sx)
    @test isapprox(test["sigma_se_dw"], sigma_se)
    @test isapprox(test["max_sx_dw"], max_sx)
    @test isapprox(test["max_sy_dw"], max_sy)
end
