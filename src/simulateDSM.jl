using ControlSystems
using LinearAlgebra

"""
    v, xn, xmax, y = simulateDSM(u, ABCD, nlev=2, x0=0)
    v, xn, xmax, y = simulateDSM(u, ntf, nlev=2, x0=0)

Compute the output of a general delta-sigma modulator with input `u`, a
structure described by `ABCD`, an initial state `x0` (default zero) and a
quantizer with a number of levels specified by `nlev`.

Alternatively, the modulator may be described by an `ntf`. The `ntf` is a
TransferFunction object. (The stf is assumed to be 1).

The structure that is simulated is the controllable canonical form
generated by ss().
"""
function simulateDSM(u, arg2, nlev=2, x0=0)
    if arg2 isa TransferFunction
        ntf = arg2
        order = length(ntf.matrix[1].z)

        ABCD = ss(-1/ntf)
        C = ABCD.C
        Sinv = orth([C' Matrix(I,order,order)]) / norm(C)
        S = inv(Sinv)
        C = C * Sinv
        if C[1] < 0
            S = -S
            Sinv = -Sinv
        end

        A = S * ABCD.A * Sinv
        B = S * ABCD.B
        B = [-B B]
        C = [1; zeros(order-1)]
        D = 1
    elseif arg2 isa Array
        if size(arg2,2) > 2 && size(arg2,2) == size(arg2,1)+1 # ABCD dimensions OK
            ABCD = arg2
            order = size(ABCD, 1) - 1
            A = ABCD[1:order, 1:order]
            B = ABCD[1:order, (order+1):(order+2)]
            C = ABCD[order+1, 1:order][:]
            D = ABCD[order+1, order+1]
        else
            error("ABCD argument does not have proper dimensions")
        end
    else
        error("arg2 is neither an ABCD matrix nor an NTF")
    end

    N = length(u)
    v = zeros(Int, N)
    y = zeros(N)
    xn = zeros(order, N)
    xmax = abs(x0)
    if x0 == 0
        x0 = zeros(order)
    end

    for i = 1:N
        y[i] = dot(C, x0) + D*u[i]
        v[i] = ds_quantize(y[i], nlev)
        x0 = A*x0 + B*[u[i]; v[i]]
        xn[:,i] = x0
        xmax = max.(abs.(x0), xmax)
    end

    return v, xn, xmax, y
end

# return an orthormal basis for the columns of A
function orth(A)
    return svd(A).U[:,1:rank(A)]
end

function ds_quantize(y, n)
    if iseven(n)
        v = 2 * floor(Int, 0.5 * y) + 1
    else
        v = 2 * floor(Int, 0.5 * (y + 1))
    end
    L = n - 1
    if v >= L
        v = L
    elseif v <= -L
        v = -L
    end
    return v
end
