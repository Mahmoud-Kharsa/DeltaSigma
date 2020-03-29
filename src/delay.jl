"""
    y = delay(x, n=1)

Delay signal `x` by `n` samples
"""
function delay(x, n=1)
    nx = length(x)
    if nx <= n
        return zeros(size(x))
    else
        return [zeros(n); x[1:nx-n]]
    end
end
