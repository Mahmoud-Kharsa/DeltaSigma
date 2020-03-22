using Test
using DeltaSigma
using MAT

@testset "simulateMS" begin

    #fixed parameters
    M = 16
    N = 2^14
    ntf = synthesizeNTF(6, 25, 1, 4)
    mtf1 = zpk([1],[0],1,1)
    mtf2a = zpk([1,1],[0.3,0.3],1,1)
    mtf2b = synthesizeNTF(2, 25, 1, 2)
    mtf4 = synthesizeNTF(4, 25*0.9, 1, 1.3)
    A0 = undbv(-3)
    A1 = undbv(-30)
    f = round(0.01*N)

    # test vectors
    n = 5 # number of tests
    A   = [  A0,   A0,    A1,    A0,   A0]
    mtf = [mtf1, mtf2a, mtf1, mtf2b, mtf4]

    # MATLAB output
    test = matread("resources/simulateMS.mat")

    for i = 1:n
        u = M * A[i] * sin.(2*pi*f/N * (0:N-1)')
        v, = simulateDSM(u, ntf, M+1)
        sv, sx, sigma_se, max_sx, max_sy = simulateMS(v, M, mtf[i])
        @test sv == test["sv"][i]
    end

end
