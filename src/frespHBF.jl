using Plots

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
        plt = plot(layout=(2,1), size=(650,500))

        F1 = evalF0(f1, z, phi)
        plot!(plt, f, abs.(F1), subplot=1, label="F1", line=(:dash, :blue))
        plot!(plt, f, phi*abs.(F2), subplot=1, label="F2", line=(:dot, :darkorange))
        plot!(plt, f, abs.(mag), subplot=1, label="HBF", line=(:darkorange))
        plot!(plt, subplot=1, xlims=(0,0.5), xticks=0:0.05:0.5, yticks=0:0.5:1)
        plot!(plt, subplot=1, title=msg)

        plot!(plt, f, dbv(mag), subplot=2, legend=false)
        plot!(plt, subplot=2, xlims=(0,0.5), xticks=0:0.05:0.5, yticks=-150:50:0)

        msg = @sprintf("pbr=%.1e", pbr)
        plot!(plt, ann=(0.0, -10, text(msg, 10, :left)))
        msg = @sprintf("sbr=%.0fdB", dbv(sbr))
        plot!(plt, ann=(0.5, dbv(sbr), text(msg, 10, :right)))

        display(plt)
        sleep(0.5) # slow down to see graphs
    end

    return mag, pbr, sbr
end
