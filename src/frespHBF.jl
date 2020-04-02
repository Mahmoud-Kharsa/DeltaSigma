using PyPlot

"""
    mag, pbr, sbr = frespHBF(f, f1, f2, phi, fp, msg)

Compute the frequency response, the passband ripple and the stopband ripple
of a Saramaki HBF. If `msg` is non-null, a plot is made. `fp` is the passband
edge. `phi` is used by designHBF.
"""
function frespHBF(f, f1, f2, phi=1, fp=0.2, msg="")
    if isempty(f)
        f = range(0, 0.5, length=1024)
    end
    if f1 isa Tuple # Presume that f1 is a (val, csd) tuple
        f1 = f1[1]
    end
    if f2 isa Tuple # Presume that f2 is a (val, csd) tuple
        f2 = f2[1]
    end

    Npts = length(f)
    w = 2 * pi * f
    z = exp.(im * w)
    cos_w = real.(z)

    n2 = length(f2)
    F2 = zeros(size(w))
    for i = 1:n2
        F2 = F2 + f2[i] * cos.(w * (2*i - 1))
    end
    F2 = F2 * 2
    mag = evalF1(f1, F2)

    passband = 1:floor(Int, 2*fp*(Npts - 1) + 1)
    stopband = Npts .+ 1 .- passband
    pbr = maximum(abs.(abs.(mag[passband]) .- 1))
    sbr = maximum(abs.(mag[stopband]))

    if !isempty(msg)
        clf()
        gcf().canvas.set_window_title("HBF Frequency Response")

        subplot(211)

        F1 = evalF0(f1, z, phi)
        plot(f, abs.(F1), linestyle="--", linewidth=1)
        plot(f, phi*abs.(F2), linestyle=":", linewidth=1)
        plot(f, abs.(mag), linewidth=1)
        legend(["F1", "F2", "HBF"])
        title(msg)
        grid()
        axis([0, 0.5, 0, 1.1])

        subplot(212)

        plot(f, dbv(mag))
        axis([0, 0.5, -150, 10])
        grid()

        msg = @sprintf("pbr=%.1e", pbr)
        text(0.0, -10, msg, verticalalignment="top")
        msg = @sprintf("sbr=%.0fdB", amp2db(sbr))
        text(0.5, dbv(sbr), msg, horizontalalignment="right", verticalalignment="bottom")

        sleep(0.5) # slow down to see graphs
    end

    return mag, pbr, sbr
end
