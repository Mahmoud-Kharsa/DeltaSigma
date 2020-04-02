using DSP

"""
    dbp(x) = 10 × log₁₀(|x|)

The dB equivalent of the power `x`
"""
function dbp(x)
    return pow2db.(abs.(x))
end
