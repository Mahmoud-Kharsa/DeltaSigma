using DSP

"""
    y = undbv(x)

Convert `x` from dB to a voltage
"""
function undbv(x)
    return db2amp.(x)
end
