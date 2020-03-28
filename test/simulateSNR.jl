using MAT

@testset "simulateSNR" begin

    # test vectors
    N = 4   # number of tests
    order = [  5,   8,  7,  7]
    osr   = [ 32,  64,  8,  8]
    H_inf = [1.5, 1.5,  2,  8]
    f0    = [  0, 1/8,  0,  0]
    nlev  = [  2,   2, 17, 17]
    amp = hcat((-120:10:-20)', [-15], (-10:0)')

    # MATLAB output
    test = matread("resources/simulateSNR.mat")

    for i = 1:N
        H = synthesizeNTF(order[i], osr[i], 1, H_inf[i], f0[i])
        snr, amp = simulateSNR(H, osr[i], amp, f0[i], nlev[i])
        @test isapprox(test["snr$i"], snr)
        @test test["amp$i"] == amp
    end

end
