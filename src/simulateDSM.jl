using ControlSystems
using LinearAlgebra

function orth(A)
    return svd(A).U[:,1:rank(A)]
end

function simulateDSM(u, arg2, nlev=[2], x0=0)
    nu = size(u, 1)
    nq = length(nlev)

    if arg2 isa TransferFunction
        ntf = arg2
        z, p, k = zpkdata(ntf)
        order = length(z[1])
        
        ABCD = ss(-1/ntf)
        C = ABCD.C
        Sinv = orth([C' Matrix(I,order,order)]) / norm(C)
        S = inv(Sinv)
        C = C*Sinv
        if C[1] < 0
            S = -S
            Sinv = -Sinv
        end
        A = S*ABCD.A*Sinv
        B = S*ABCD.B
        B = [-B B]
        C = [1 zeros(1,order-1)]
        D = 1
    elseif arg2 isa Array
        if size(arg2,2) > 2 && size(arg2,2) == (nu+size(arg2,1))    # ABCD dimensions OK
            ABCD = arg2
            order = size(ABCD, 1) - nq

            A = ABCD[1:order, 1:order]
            B = ABCD[1:order, (order+1):(order+nu+nq)]
            C = ABCD[(order+1):(order+nq), 1:order]
            D = ABCD[(order+1):(order+nq), (order+1):(order+nu)]
        else
            throw(ArgumentError("ABCD argument does not have proper dimensions"))
        end
    else
        throw(ArgumentError("second argument is neither an ABCD matrix nor an NTF"))
    end

    N = length(u)
    v = zeros(nq, N)
    y = zeros(nq, N)
    xn = zeros(order, N)
    xmax = abs(x0)
    if x0 == 0
        x0 = zeros(order, 1)
    end

    for i = 1:N
        y[:,i] = C*x0 + D*u[:,i]
        v[:,i] = ds_quantize(y[:,i], nlev)
        x0 = A*x0 + B*[u[:,i]; v[:,i]]
        xn[:,i] = x0
        xmax = max.(abs.(x0), xmax)
    end

    return v, xn, xmax, y
end

function ds_quantize(y, n)
    v = zeros(size(y))

    i = iseven.(n)
    v[i,:] .= 2 * floor.(0.5 * y[i,:]) .+ 1

    i = isodd.(n)
    v[i,:] .= 2 * floor.(0.5 * (y[i,:] .+ 1))

    for qi = 1:length(n)
        L = n[qi] - 1

        i = v[qi,:] .> L
        v[qi,i] .= L

        i = v[qi,:] .< -L
        v[qi,i] .= -L
    end

    return v
end
