using DSP
using FFTW
using MAT

@testset "calculateSNR" begin
    # MATLAB output
    test = matread("resources/calculateSNR.mat")

    # low pass
    u = 0.5 * sin.(2*pi*85/8192 * (0:8191))
    ntf = synthesizeNTF(5, 32, 1)
    v, = simulateDSM(u, ntf)
    spec = fft(v .* hanning(8193)[1:8192]/2048)
    snr = calculateSNR(spec[3:129], 83)
    @test isapprox(test["snr_low"], snr)

    # band pass
    ntf = synthesizeNTF(8, 64, 1, 1.5, 0.125)
    u = 0.5 * sin.(2*pi*1045/8192 * (0:8191))
    v, = simulateDSM(u, ntf)
    spec = fft(v .* hanning(8193)[1:8192]/2048)
    snr = calculateSNR(spec[992:1056], 54)
    @test isapprox(test["snr_band"], snr)
end
