using LinearAlgebra

"""
    g = rmsGain(H, f1, f2, N=100)

Compute the root mean-square gain of the discrete-time
tf `H` in the frequency band (`f1`, `f2`)
"""
function rmsGain(H, f1, f2, N=100)
    w = range(2*pi*f1, 2*pi*f2, length=N)
    return norm(evalTF(H, exp.(im*w))) / sqrt(N)
end
