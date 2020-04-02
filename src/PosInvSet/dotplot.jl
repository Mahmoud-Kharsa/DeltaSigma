using PyPlot

"""
    dotplot(p, color, marker)

Plot dots with `color` and `marker`
"""
function dotplot(p, color, marker)
    if marker == "o" || marker == "s"
        plot(p[1,:], p[2,:], color=color, linestyle="None", marker=marker, markerfacecolor="none", markeredgewidth=0.5)
    else
        plot(p[1,:], p[2,:], color=color, linestyle="None", marker=marker, ms=2)
    end
end
