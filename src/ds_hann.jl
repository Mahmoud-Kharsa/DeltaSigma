"""
    w = ds_hann(n)

A Hann window of length `n`. Does not smear tones located exactly in a bin.
"""
function ds_hann(n)
    return 0.5 * (1 .- cos.(2*pi*(0:n-1)/n))
end
