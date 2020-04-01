"""
    y = bquantize(x, nsd=3, abstol=eps, reltol=10*eps)

Bidirectionally quantize a n by 1 vector x to nsd signed digits, terminate early
if the error is less than the specified tolerances.

`y` is a 2-tuple of arrays with the following fields:
- `y`[1][i] is the quantized value in floating-point form
- `y`[2][i] is a 2-by-nsd (or less) matrix containing the powers of two
    (first row) and their signs (second row).

See also bunquantize.
"""
function bquantize(x, nsd=3, abstol=eps(Float64), reltol=10*eps(Float64))
    n = length(x)
    offset = -log2(0.75)

    y_val = zeros(n)
    y_csd = Array{Matrix{Int}}(undef, n)

    for i = 1:n
        xp = x[i]
        y_val[i] = 0.0
        for j = 1:nsd
            error = abs(y_val[i] - x[i])
            if error <= abstol || error <= abs(x[i])*reltol
                break
            end
            p = floor(Int, log2(abs(xp)) + offset)
            p2 = 2.0 .^ p
            sx = cmp(xp, 0)
            xp = xp - sx*p2
            y_val[i] = y_val[i] + sx*p2
            if isassigned(y_csd, i)
                y_csd[i] = hcat(y_csd[i], [p; sx])
            else
                y_csd[i] = reshape([p; sx], (2,1))
            end
        end
    end

    return y_val, y_csd
end
