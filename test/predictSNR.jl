using MAT

@testset "predictSNR" begin
    # MATLAB output
    test = matread("resources/predictSNR.mat")

    # low pass
    ntf = synthesizeNTF(5, 32, 1)
    snr, amp, k0, k1, sigma_e2 = predictSNR(ntf, 32)
    @test isapprox(vec(test["snr_low"]), snr)
    @test [-120:10:-20; -15; -10:0] == amp   # default
    @test isapprox(vec(test["k0_low"]), k0)
    @test isapprox(vec(test["k1_low"]), k1)
    @test isapprox(vec(test["sigma_e2_low"]), sigma_e2)

    # band pass
    ntf = synthesizeNTF(8, 64, 1, 1.5, 0.125)
    amp_test = [-120:10:-20; -15; -10:0]
    snr, amp, k0, k1, sigma_e2 = predictSNR(ntf, 64, amp_test, 0.125)
    # Allow higher error tolerance because of differences in MATLAB and Julia
    # interpolation libraries. Default tolerance is sqrt(eps).
    rtol = 1e4*sqrt(eps(Float64))
    @test isapprox(vec(test["snr_band"]), snr, rtol=rtol)
    @test amp_test == amp   # passed through without changes
    @test isapprox(vec(test["k0_band"]), k0, rtol=rtol)
    @test isapprox(vec(test["k1_band"]), k1, rtol=rtol)
    @test isapprox(vec(test["sigma_e2_band"]), sigma_e2, rtol=rtol)
end
