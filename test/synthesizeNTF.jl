using Test
using ControlSystems
using DeltaSigma
using MAT

@testset "synthesizeNTF" begin

    # error tolerance
    rtol = 10 * eps(Float64)

    # test vectors
    N = 10   # number of tests
    order = [  5,   5,     8,     8, 7, 7,   5,  6,  2,    4]
    osr   = [ 32,  32,    64,    64, 8, 8,  42, 25, 25, 22.5]
    opt   = [  0,   1,     2,     1, 1, 1,   1,  1,  1,    1]
    H_inf = [1.5, 1.5,   1.5,   1.5, 2, 8, 1.5,  4,  2,  1.3]
    f0    = [  0,   0, 0.125, 0.125, 0, 0,   0,  0,  0,    0]

    # MATLAB output
    test = matread("synthesizeNTF.mat")

    for i = 1:N
        ntf = synthesizeNTF(order[i], osr[i], opt[i], H_inf[i], f0[i])
        data = zpkdata(ntf)
        z, p, k = data[1][1], data[2][1], data[3][1]
        testz = vec(test["ntf$i"]["z"])
        testp = vec(test["ntf$i"]["p"])
        testk = vec(test["ntf$i"]["k"])
        sort!(z, by = x->(real(x), imag(x)))
        sort!(p, by = x->(real(x), imag(x)))
        sort!(testz, by = x->(real(x), imag(x)))
        sort!(testp, by = x->(real(x), imag(x)))
        @test isapprox(testz, z; rtol=rtol)
        @test isapprox(testp, p; rtol=rtol)
        @test isapprox(testk, k; rtol=rtol)
    end

end