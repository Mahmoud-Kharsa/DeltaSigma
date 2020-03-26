"""
    A, B, C, D = partitionABCD(ABCD, m)

Partition `ABCD` into `A`, `B`, `C`, `D` for an `m`-input state-space system.
"""
function partitionABCD(ABCD, m=NaN)
    if isnan(m)
        n = minimum(size(ABCD)) - 1
        m = size(ABCD, 2) - n
    else
        n = size(ABCD, 2) - m
    end
    r = size(ABCD, 1) - n

    A = ABCD[1:n, 1:n]
    B = ABCD[1:n, n+1:n+m]
    C = ABCD[n+1:n+r, 1:n]
    D = ABCD[n+1:n+r, n+1:n+m]

    return A, B, C, D
end
