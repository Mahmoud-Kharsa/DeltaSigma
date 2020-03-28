using DSP
using FFTW

"""
    snr, amp = simulateSNR(ntf, osr, amp, f0=0, nlev=2, f=1/(4*osr), k=13)
    snr, amp = simulateSNR(ABCD, osr, amp, f0=0, nlev=2, f=1/(4*osr), k=13)

Determine the `snr` for a delta-sigma modulator by using simulations.

The modulator is described by a noise transfer function `ntf` and the
number of quantizer levels `nlev`. Alternatively, the first argument may
be an `ABCD` matrix.

The band of interest is defined by the oversampling ratio `osr` and the
center frequency `f0`.

The input signal is characterized by the `amp` vector and the `f` variable.
- `amp` defaults to [-120, -110, ... -20, -15, -10, -9, -8, ... 0] dB, where
  0 dB means a full-scale (peak value = `nlev`-1) sine wave.
- `f` is the input frequency, normalized such that 1 -> fs; `f` is rounded
  to an FFT bin.

Using sine waves located in FFT bins, the `snr` is calculated as the ratio
of the sine wave power to the power in all in-band bins other than those
associated with the input tone. Due to spectral smearing, the input tone
is not allowed to lie in bins 0 or 1. The length of the FFT is 2`ᵏ`.
"""
function simulateSNR(arg1, osr, amp=[-120:10:-20; -15; -10:0], f0=0, nlev=2, f=NaN, k=13)
    osr_mult = 2
    if f0 != 0
        osr_mult = 4
    end
    if isnan(f)
        f = f0 + 0.5/(osr*osr_mult)   # Halfway across the band
    end
    M = nlev - 1

    if abs(f-f0) > 1/(osr*osr_mult)
        @warn "the input tone is out-of-band."
    end

    N = 2^k
    if N < 8*2*osr  # Require at least 8 bins to be "in-band"
        @warn "Increasing k to accommodate a large oversampling ratio."
        k = ceil(Int, log2(8*2*osr))
        N = 2^k
    end

    F = round(Int, f*N)
    if abs(F) <= 1
        @warn "Increasing k to accomodate a low input frequency."
        # Want f*N >= 1
        k = ceil(Int, log2(1/f))
        N = 2^k
        F = 2
    end

    Ntransient = 100
    soft_start = 0.5 * (1 .- cos.(2*pi/Ntransient * (0:(Ntransient/2-1))))
    tone = M * sin.(2*pi*F/N * (0:(N+Ntransient-1)))
    tone[1:(Ntransient÷2)] = tone[1:(Ntransient÷2)] .* soft_start
    window = hanning(N+1)[1:N]  # Hann window

    if f0 == 0
        # Exclude DC and its adjacent bin
        inBandBins = N÷2 .+ (3:round(Int, N/(osr_mult*osr)))
        F = F - 2
    else
        f1 = round(Int, N * (f0 - 1/(osr_mult*osr)))
        inBandBins = N÷2 .+ (f1:round(Int, N * (f0 + 1/(osr_mult*osr)))) # Should exclude DC
        F = F - f1 + 1
    end

    i = 1
    snr = zeros(size(amp))
    for A = 10 .^ (amp/20)
        v, = simulateDSM(A*tone, arg1, nlev)
        hwfft = fftshift(fft(window .* v[(1+Ntransient):(N+Ntransient)]))
        snr[i] = calculateSNR(hwfft[inBandBins], F)
        i = i+1
    end

    return snr, amp
end
