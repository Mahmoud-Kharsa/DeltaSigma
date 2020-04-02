using Polynomials

"""
    y = evalRPoly(roots, x, k=1)

Compute the value of a polynomial which is given in terms of its `roots`.
"""
function evalRPoly(roots, x, k=1)
    f = poly(roots) * k
    return f(x)
end
