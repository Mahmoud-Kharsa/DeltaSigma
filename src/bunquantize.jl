"""
    y = bunquantize(q)

Calculate the value corresponding to a bidirectionally quantized quantity.
`q` is a 2n by m matrix containing the powers of 2 and their signs for each
quantized value.
"""
function bunquantize(q)
    n = size(q, 1) ÷ 2
    y = zeros(n)
    signs = 2:2:2*n
    powers = signs .- 1
    for i = 1:size(q, 2)
        y = y + 2.0 .^ q[powers,i] .* q[signs,i]
    end
    return y
end
