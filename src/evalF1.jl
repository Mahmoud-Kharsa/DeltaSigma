"""
    F1 = evalF1(f1, z, phi)

Calculate the values of the `F1` filter (tranformed prototype filter) of a
Saramaki HBF at the given points.
"""
function evalF1(f1, z, phi=1.0)
    z = z / phi
    F1 = 0.5 * ones(length(z))
    for i = 1:length(f1)
        F1 = F1 + f1[i] * z .^ (2*i - 1)
    end
    return F1
end
