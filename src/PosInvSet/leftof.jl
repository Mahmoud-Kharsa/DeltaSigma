"""
    y = leftof(p, a, b)

Return true if the point `p` is to the left of the line `a``b`.
For a n×2 list of points `p`, return a vector of results.
"""
function leftof(p, a, b)
    # Translate to the origin and do the check
    p = p - ones(size(p, 1)) * a
    b = b - a
    return b[1]*p[:,2] .> b[2]*p[:,1]
end
