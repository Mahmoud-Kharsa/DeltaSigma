"""
    h = evalTF(tf, z)

Evaluates the transfer function described by `tf` at the point(s) given in the
`z` vector. `tf` must be a TransferFunction object.
"""
function evalTF(tf, z)
    if size(z) == ()
        return tf(z)[1, 1]
    else
        return tf(z, false)[:, 1, 1]
    end
end
