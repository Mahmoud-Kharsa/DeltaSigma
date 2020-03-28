using MAT

@testset "peakSNR" begin
    # MATLAB output
    test = matread("resources/peakSNR.mat")

    # test vectors
    osr  = [32,    64,  8,  8]
    f0   = [ 0, 0.125,  0,  0]
    nlev = [ 2,     2, 17, 17]
    ntf = [
        synthesizeNTF(5, osr[1], 1, 1.5, f0[1]),
        synthesizeNTF(8, osr[2], 1, 1.5, f0[2]),
        synthesizeNTF(7, osr[3], 1,   2, f0[3]),
        synthesizeNTF(7, osr[4], 1,   8, f0[4])
    ]
    amp_test = [-120:10:-20; -15; -10:0]

    for i = 1:4
        snr, amp = simulateSNR(ntf[i], osr[i], amp_test, f0[i], nlev[i])
        peak_snr, peak_amp = peakSNR(snr, amp)
        @test isapprox(test["peak_snr_$i"], peak_snr)
        @test test["peak_amp_$i"] == peak_amp
    end
end
