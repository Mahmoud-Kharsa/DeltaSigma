using Statistics

"""
    peak_snr, peak_amp = peakSNR(snr, amp)

Estimate the snr peak by threading a line of slope 1 through the (`amp`,`snr`)
data. Both `amp` and `snr` are expressed in dB.
"""
function peakSNR(snr, amp)
    # Delete garbage data
    i = abs.(snr) .!= Inf
    snr = snr[i]
    amp = amp[i]
    i = .!isnan.(snr)
    snr = snr[i]
    amp = amp[i]
    i = snr .>= 3
    snr = snr[i]
    amp = amp[i]

    n = length(amp)

    # Sort by amplitude
    i = sortperm(amp)
    amp = amp[i]
    snr = snr[i]

    i = 1
    m = 0
    while any(i .!= 0) && n > 3
        # Draw a 45-degree median line through the data
        tmp = sort(snr-amp)
        if rem(n,2) == 0
            m = mean(tmp[n÷2 .+ [1,0]])
        else
            m = tmp[(n+1)÷2]
        end
        # Discard data that is more than 6dB off
        i = abs.(amp - snr .+ m) .> 6
        snr = snr[.!i]
        amp = amp[.!i]
        n = length(amp)
    end

    peak_amp = maximum(amp)
    peak_snr = peak_amp + m
    return peak_snr, peak_amp
end
