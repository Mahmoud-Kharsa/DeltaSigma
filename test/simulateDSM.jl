using MAT

@testset "simulateDSM" begin
    # MATLAB output
    test = matread("resources/simulateDSM.mat")

    # ntf input test
    u1 = 0.5 * sin.(2*pi*85/8192 * (0:8191))
    u2 = 0.5 * sin.(2*pi*1045/8192 * (0:8191))
    u3 = 8 * sin.(2*pi*146/8192 * (0:8191))
    u4 = 16 * undbv(-3) * sin.(2*pi*164/16384 * (0:16383))
    u5 = 16 * undbv(-30) * sin.(2*pi*164/16384 * (0:16383))
    ntf1 = synthesizeNTF(5, 32, 1)
    ntf2 = synthesizeNTF(8, 64, 1, 1.5, 0.125)
    ntf3 = synthesizeNTF(7, 8, 1, 2)
    ntf4 = synthesizeNTF(7, 8, 1, 8)
    ntf5 = synthesizeNTF(6, 25, 1, 4)
    u    = [  u1,   u2,   u3,   u3,   u4,   u5]
    ntf  = [ntf1, ntf2, ntf3, ntf4, ntf5, ntf5]
    nlev = [   2,    2,   17,   17,   17,   17]
    for i = 1:length(u)
        v, _, _, y = simulateDSM(u[i], ntf[i], nlev[i])
        @test vec(test["v_$i"]) == v
        # Allow higher error tolerance because of differences in MATLAB and
        # Julia state space models. This also why xn and xmax are ignored;
        # the x vectors aren't the same. Default tolerance is sqrt(eps).
        @test isapprox(vec(test["y_$i"]), y, rtol=10*sqrt(eps(Float64)))
    end

    # ABCD input test
    u = range(0, 0.6; length=30)
    ntf = synthesizeNTF(5, 42, 1)
    a, g, b, c = realizeNTF(ntf)
    ABCD = stuffABCD(a, g, b, c)
    for i = 1:length(u)
        v, xn, xmax, y = simulateDSM(u[i]*ones(10000), ABCD)
        @test vec(test["v"][i:i,:]) == v
        @test isapprox(test["xn"][i,:,:], xn)
        @test isapprox(test["xmax"][i,:,:], xmax)
        @test isapprox(test["y"][i,:,:], y)
    end
end
