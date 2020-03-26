using ControlSystems
using LinearAlgebra

"""
    ntf, stf = calculateTF(ABCD, k=1)

Calculate the NTF and STF of a delta-sigma modulator whose loop filter
is described by the `ABCD` matrix, assuming a quantizer gain of `k`.

The `ntf` and `stf` are TransferFunction objects.
"""
function calculateTF(ABCD, k=[1.0])
    nq = length(k)
    n = size(ABCD, 1) - nq
    nu = size(ABCD, 2) - size(ABCD, 1)
    A, B, C, D = partitionABCD(ABCD, nu+nq)
    B1 = B[:,1:nu]
    B2 = B[:,nu+1:end]
    D1 = D[:,1:nu]
    if any(x -> x != 0, D[:,nu+1:end])
        error("D must be zero")
    end
    K = diagm(k)

    # Find the noise transfer function by forming the closed-loop
    # system (sys_cl) in state-space form.
    Acl = A + B2*K*C
    Bcl = [(B1 + B2*K*D1) B2]
    Ccl = K*C
    Dcl = [K*D1 diagm(ones(nq))]

    sys_cl = ss(Acl, Bcl, Ccl, Dcl, 1)
    tfs = zpk(sys_cl)
    tol = min(1e-3, max(1e-6, eps(Float64)^(1/(size(ABCD,1)))))
    mtfs = minreal(tfs, tol)
    stf = TransferFunction(mtfs.matrix[:,1:nu], mtfs.Ts)
    ntf = TransferFunction(mtfs.matrix[:,nu+1:nu+nq], mtfs.Ts)

    return ntf, stf
end
