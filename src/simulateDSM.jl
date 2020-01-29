using ControlSystems
using LinearAlgebra

function orth(A)
    return svd(A).U[:,1:rank(A)]
end

function simulateDSM(u, ntf, nlev=[2], x0=0)
    nu = size(u, 1)
    nq = length(nlev)

    z, p, k = zpkdata(ntf)
    order = length(z[1])
    
    ABCD = ss(-1/ntf)
    A = ABCD.A
    B2 = ABCD.B
    C = ABCD.C
    D2 = ABCD.D

    Sinv = orth([C' Matrix(I,order,order)]) / norm(C)
    S = inv(Sinv)
    C = C*Sinv
    if C[1] < 0
        S = -S
        Sinv = -Sinv
    end
    A = S*A*Sinv
    B2 = S*B2
    C = [1 zeros(1,order-1)]
    D2 = 0
    B1 = -B2
    D1 = 1
    B = [B1 B2]

    N = length(u)
    v = zeros(nq, N)
    y = zeros(nq, N)
    xn = zeros(order, N)
    xmax = abs(x0)
    if x0 == 0
        x0 = zeros(order, 1)
    end

    for i = 1:N
        y[:,i] = C*x0 + D1*u[:,i]
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
