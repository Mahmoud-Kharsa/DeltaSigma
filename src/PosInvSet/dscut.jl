"""
    p = dscut(p1, y1, p2, y2)

return the point `p` at which y = 0, assuming y varies linearly
from `y1` at `p1` to `y2` at `p2`.
"""
function dscut(p1, y1, p2, y2)
    k1 = y2 / (y2 - y1)
    k2 = 1 - k1
    n = size(p1, 1)
    return k1 .* p1 + k2 .* p2
end
