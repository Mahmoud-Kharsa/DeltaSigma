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
    dbg = (length(options) >= 1) ? options[1] : 0
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

        if dbg == 1
            plot(legend=false, size=(550,450))
            plot!(xlims=(axis1[1], axis1[2]), ylims=(axis1[3], axis1[4]))
            plot!(xticks=round(axis1[1]):0.5:round(axis1[2]))
            plot!(yticks=round(axis1[3]):1:round(axis1[4]))

            plot!(ec[1,:], ec[2,:], line=(:scatter), marker=(6, :white, stroke(:red)))
            plot!( x[1,:], x[2,:], line=(:scatter), marker=(2, :black))

            outi = (out .!= 0)[:]
            plot!(ns[1,outi], ns[2,outi], line=(:scatter), marker=(:square, :white, stroke(:red)))

            polyplot!(s, :blue)
            polyplot!(s1, :magenta)
            polyplot!(s2, :cyan)

            str = @sprintf("Iteration %d: %d image vertices outside", i, sum(outi))
            plot!(title=str)
            display(plot!())
            sleep(0.1)
        end

        if all(out .== 0)
            break
        end

        s = hull2d(ns')[1]'
    end

    return s
end
