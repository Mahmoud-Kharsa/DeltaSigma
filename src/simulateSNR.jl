function simulateSNR(arg1, osr; amp=vcat(-120:10:-20, -15, -10:0), f0=0, nlev=2, f=NaN, k=13)
    osr_mult = 2
    if f0 != 0
        osr_mult = 4
    end
    if isnan(f)
        f = f0 + 0.5/(osr*osr_mult)   # Halfway across the band
    end
    M = nlev - 1

    if abs(f-f0) > 1/(osr*osr_mult)
        println(stderr, "Warning: the input tone is out-of-band.")
    end

    N = 2^k
    if N < 8*2*osr  # Require at least 8 bins to be "in-band"
        println("Warning: Increasing k to accommodate a large oversampling ratio.")
        k = ceil(Int, log2(8*2*osr))
        N = 2^k
    end
    F = round(Int, f*N)
    if abs(F) <= 1
        println("Warning: Increasing k to accomodate a low input frequency.")
        # Want f*N >= 1
        k = ceil(Int, log2(1/f))
        N = 2^k
        F = 2
    end

    Ntransient = 100
    soft_start = 0.5 * (1 .- cos.(2*pi/Ntransient * (0:(Ntransient/2-1))))
    tone = M * sin.(2*pi*F/N * (0:(N+Ntransient-1)))
    tone[1:div(Ntransient,2)] = tone[1:div(Ntransient,2)] .* soft_start
    window = 0.5 * (1 .- cos.(2*pi*(0:N-1)/N))  # Hann window
    if f0 == 0
        # Exclude DC and its adjacent bin
        inBandBins = div(N,2) .+ (3:round(Int, N/(osr_mult*osr)))
        F = F-2
    else
        f1 = round(Int, N * (f0 - 1/(osr_mult*osr)))
        inBandBins = div(N,2).+ (f1:round(Int, N * (f0 + 1/(osr_mult*osr)))) # Should exclude DC
        F = F-f1+1;
    end
    
    snr = zeros(size(amp))
    i = 1
    for A = 10 .^ (amp/20)
        v, = simulateDSM(A*tone, arg1, nlev)
        hwfft = fftshift(fft(window .* v[(1+Ntransient):(N+Ntransient)]))
        snr[i] = calculateSNR(hwfft[inBandBins], F)
        i = i+1;
    end

    return snr, amp
end