using Test
using DeltaSigma
using MAT

@testset "simulateDSM" begin

    # MATLAB output
    test = matread("resources/simulateDSM.mat")

    # transfer function input test
    N = 6   # number of tests
    A     = [ 0.5,  0.5,    8,    8, 16*undbv(-3), 16*undbv(-30)]
    f     = [  85, 1045,  146,  146,          164,           164]
    n     = [8192, 8192, 8192, 8192,        16384,         16384]
    order = [   5,    8,    7,    7,            6,             6]
    osr   = [  32,   64,    8,    8,           25,            25]
    H_inf = [ 1.5,  1.5,    2,    8,            4,             4]
    f0    = [   0,  1/8,    0,    0,            0,             0]
    nlev  = [   2,    2,   17,   17,           17,            17]
    for i = 1:N
        u = A[i] * sin.(2*pi*f[i]/n[i] * (0:n[i]-1)')
        H = synthesizeNTF(order[i], osr[i], 1, H_inf[i], f0[i])
        v, = simulateDSM(u, H, nlev[i])
        @test test["v$i"] == v
    end

    # ABCD input test
    u = range(0, 0.6; length=30)
    ABCD = test["ABCD"]  # TODO: generate using realizeNTF and stuffABCD
    for i = 1:30
        v, xn, xmax, = simulateDSM(u[i]*ones(1,10000), ABCD)
        @test test["v"][i:i,:] == v
        @test isapprox(test["xn"][i,:,:], xn)
        @test isapprox(test["xmax"][i,:,:], xmax)
    end

end
