using LaTeXStrings
using PyPlot
using Printf

function dsdemo7()
    println("Invariant Set for MOD2 (find2dPIS)")

    A = [1 0; 1 1]
    B = [1 -1; 1 -2]
    C = [0 1]
    D = [0 0]
    ABCD = [A B; C D]
    order = 2

    fig = figure(1)
    fig.clf()
    fig.canvas.set_window_title("Working Invariant Set")

    u = 1 / pi
    s = find2dPIS(u, ABCD, true)

    N = 10000
    skip = 100
    _, x, = simulateDSM(u * ones(N+skip), ABCD, 2)
    x = x[:,1+skip:N+skip]
    nv = size(s, 2)
    splus, eplus, sminus, eminus = dssplit2d(u, ABCD, s)
    Buv = B * [u, 1]
    s1 = A*splus + Buv*ones(1,size(splus, 2))
    Buv = B * [u, -1]
    s2 = A*sminus + Buv*ones(1,size(sminus, 2))
    ns = [s1 s2]
    out = outconvex2d(ns, s)

    fig = figure(2)
    fig.clf()
    fig.canvas.set_window_title("Final Invariant Set")

    grid()
    dotplot(x, "k", ".")
    polyplot(s, "b")
    polyplot(s1, "m")
    polyplot(s2, "c")
    outi = (out .!= 0)[:]
    dotplot(ns[:,outi], "r", "s")
    str = @sprintf("Final Object: %d image vertices outside", sum(outi))
    title(str)
    xlabel(L"x_1")
    ylabel(L"y_1")

    @printf("%d points from the %d simulated states are outside.\n", sum(outconvex2d(x, s)), N)
    @printf("%d image points are outside.\n", sum(out))
    @printf("The returned polygon has %d vertices.\n", size(s, 2))

    return nothing
end
