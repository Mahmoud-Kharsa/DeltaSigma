using DSP

"""
    y = undbp(x)

Convert `x` from dB to a power
"""
function undbp(x)
    return db2pow.(x)
end
