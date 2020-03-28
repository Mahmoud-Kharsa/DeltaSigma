"""
    v, i = hull2d(p)

Convex hull of a set of 2D points.
`p` is a (np)×2 array; `v` is a (nv+1)×2 array, with the last point being a
duplicate of the first, (that way, plot(`v`[:,1], `v`[:,2]) yields a plot of
the polygon) and `i` is the index of each vertex in the `p` array. `v` starts
at the lower leftmost point and moves ccw around the hull.

The algorithm was created by Richard Schreier in the MATLAB Delta-Sigma Toolbox.
"""
function hull2d(p)
    # Find the (lower)left-most points.
    x = minimum(p[:,1])
    leftmost = findall(a -> a == x, p[:,1])
    y, il = findmin(p[leftmost,2])
    il = leftmost[il]
    l = [x y]

    # Find the (upper)right-most points.
    x = maximum(p[:,1])
    rightmost = findall(a -> a == x, p[:,1])
    y, ir = findmax(p[rightmost,2])
    ir = rightmost[ir]
    r = [x y]

    if l == r		# Degenerate cases
        return l, il
    end

    # Split the points into those above and those below the l-r line.
    isAbove = leftof(p, l, r)
    ia = findall(isAbove)
    ib = findall(.!isAbove)
    above = p[ia,:]
    below = p[ib,:]

    # Sort them in terms of increasing first coordinate.
    if any(isAbove)
        isort = sortperm(above[:,1])
        above = [l; above[isort, :]]	# l must be first.
        ia = [il; ia[isort]]
    else    # no points above
        above = l
        ia = il
    end
    isort = sortperm(below[:,1])
    below = [below[isort,:]; r] # includes the l and r points r must be last.
    ib = [ib[isort]; ir]

    # Move along the underside, building the vertex list as we go.
    a = below[1,:]'
    nb = size(below, 1)
    b = below[2,:]'
    v = [a; b]
    i = [ib[1]; ib[2]]
    nv = 2

    for n = 3:nb
        p = below[n,:]'
        while !leftof(p, a, b)[1]   # backtrack vertices
            nv = nv - 1
            v = v[1:nv,:]
            i = i[1:nv]
            b = a
            if nv > 1
                a = v[nv-1,:]'
            else
                break
            end
        end
        v = [v; p]
        i = [i; ib[n]]
        nv = nv + 1
        a = b
        b = p
    end

    # Move along the top side, continuing to build the vertex list.
    na = size(above, 1)
    for n = na:-1:1
        p = above[n,:]'
        while !leftof(p, a, b)[1] && nv > 2
            nv = nv - 1
            v = v[1:nv,:]
            i = i[1:nv]
            b = a
            a = v[nv-1,:]'
        end
        v = [v; p]
        i = [i; ia[n]]
        nv = nv + 1
        a = b
        b = p
    end

    return v, i
end
