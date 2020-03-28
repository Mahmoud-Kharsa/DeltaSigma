using ControlSystems
using MAT

@testset "synthesizeNTF" begin
    # MATLAB output
    test = matread("resources/synthesizeNTF.mat")

    # test vectors
    order = [  5,   5,     8,     8, 7, 7,   5,  6,  2,    4]
    osr   = [ 32,  32,    64,    64, 8, 8,  42, 25, 25, 22.5]
    opt   = [  0,   1,     2,     1, 1, 1,   1,  1,  1,    1]
    H_inf = [1.5, 1.5,   1.5,   1.5, 2, 8, 1.5,  4,  2,  1.3]
    f0    = [  0,   0, 0.125, 0.125, 0, 0,   0,  0,  0,    0]

    for i = 1:length(order)
        ntf = synthesizeNTF(order[i], osr[i], opt[i], H_inf[i], f0[i])
        data = zpkdata(ntf)

        z = ntf.matrix[1].z
        p = ntf.matrix[1].p
        k = ntf.matrix[1].k
        sort!(z, by=x->(real(x), imag(x)))
        sort!(p, by=x->(real(x), imag(x)))

        testz = vec(test["ntf_$i"]["z"])
        testp = vec(test["ntf_$i"]["p"])
        testk = test["ntf_$i"]["k"]
        sort!(testz, by=x->(real(x), imag(x)))
        sort!(testp, by=x->(real(x), imag(x)))

        @test isapprox(testz, z)
        @test isapprox(testp, p)
        @test isapprox(testk, k)
    end
end
