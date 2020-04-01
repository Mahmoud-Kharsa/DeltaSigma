"""
    F0 = evalF0(f1, z, phi)

Calculate the values of the `F0` (prototype) filter of a Saramaki HBF
at the given points.
"""
function evalF0(f1, z, phi)
    return evalF1(f1, 0.5*(z + 1 ./ z), phi)
end
