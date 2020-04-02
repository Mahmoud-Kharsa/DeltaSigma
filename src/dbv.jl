using DSP

"""
    dbv(x) = 20 × log₁₀(|x|)

The dB equivalent of the voltage `x`
"""
function dbv(x)
    return amp2db.(abs.(x))
end
