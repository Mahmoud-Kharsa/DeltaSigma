using LinearAlgebra

"""
    out = outconvex2d(x, p)

Test if each of the `x` points are inside the convex polygon `p`, and return the
number of inequalities failed by each point. `p` is a 2×n counter-clockwise list
of vertices, with the first vertex duplicated.
"""
function outconvex2d(x, p)
    n = size(p, 2)

    # form A,B such that internal points satisfy Ax <= B
    A = [p[2,2:n] - p[2,1:n-1] p[1,1:n-1] - p[1,2:n]];
    B = zeros(n-1)
    for i = 1:n-1
        B[i] = dot(A[i,:], p[:,i])
    end

    return sum(A*x .> B*ones(1, size(x,2)), dims=1)
end
