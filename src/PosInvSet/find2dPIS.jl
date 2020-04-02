using PyPlot

"""
    s = find2dPIS(u, ABCD, options)

Find a positively invariant set for the 2nd-order binary modulator whose
loop filter is described by `ABCD` and whose input is a constant `u`.
Return `s` = [] if no invariant set is found.

`options` = [dbg, itnLimit, expFactor, N, skip] defaults to [0, 100, 0.005, 1000, 100],
dbg=1 causes debugging information to be displayed.
"""
function find2dPIS(u, ABCD, options=[])
    n = size(ABCD, 1) - 1
    if n != 2
        error("The modulator must be second order.")
    end

    nlev = 2
    dbg = (length(options) >= 1) ? options[1] : false
    itnLimit = (length(options) >= 2) ? options[2] : 100
    expFactor = (length(options) >= 3) ? options[3] : 0.005
    N = (length(options) >= 4) ? options[4] : 1000
    skip = (length(options) >= 5) ? options[5] : 100

    # Compute a few iterations of difference equations
    if size(u) == ()
        un = u * ones(N+skip)
    else
        error("Argument 1 (u) has the wrong dimensions.")
    end
    v, x, = simulateDSM(un, ABCD, nlev)
    x = x[:,skip+1:skip+N]

    xmin = minimum(x[1,:])
    xmax = maximum(x[1,:])
    dx = xmax - xmin

    ymin = minimum(x[2,:])
    ymax = maximum(x[2,:])
    dy = ymax - ymin

    axis1 = [xmin-dx/4, xmax+dx/4, ymin-dy/4, ymax+dy/4]

    # Take the convex hull of the result
    s = hull2d(x')[1]'
    ec = mean(x, dims=2)

    for i = 1:itnLimit
        # Inflate the hull
        shift = ec * ones(1, size(s, 2))
        s = shift + (1 + expFactor) * (s - shift)

        # Split the set
        splus, eplus, sminus, eminus = dssplit2d(u, ABCD, s)
        # Map the two halves
        s1 = dsmap(u, ABCD, 2, splus, eplus, 1)
        s2 = dsmap(u, ABCD, 2, sminus, eminus, -1)
        ns = [s1[:,1:size(s1,2)-1] s2[:,1:size(s2,2)-1]]

        # Test for inclusion: ns inside s (the inflated hull)
        out = outconvex2d(ns, s)

        if dbg
            clf()
            grid()
            axis(axis1)
            dotplot(x, "k", ".")
            dotplot(ec, "r", "o")
            polyplot(s, "b")
            polyplot(s1, "m")
            polyplot(s2, "c")
            outi = (out .!= 0)[:]
            dotplot(ns[:,outi], "r", "s")
            str = @sprintf("Iteration %d: %d image vertices outside", i, sum(outi))
            title(str)
            sleep(0.1)
        end

        if all(out .== 0)
            break
        end

        s = hull2d(ns')[1]'
    end

    return s
end
