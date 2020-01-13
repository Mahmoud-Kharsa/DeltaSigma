function calculateSNR(hwfft, f, nsig=1)
    signalBins = f-nsig+1:f+nsig+1
    signalBins = signalBins[signalBins .> 0]
    signalBins = signalBins[signalBins .<= length(hwfft)]
    s = norm(hwfft[signalBins])

    noiseBins = collect(1:length(hwfft))
    deleteat!(noiseBins, signalBins)
    n = norm(hwfft[noiseBins])

    if n == 0
        snr = Inf
    else
        snr = dbv(s/n)
    end

    return snr
end
