using Plots

"""
    polyplot!(p, color)

Plot a polygon given by the point list `p` with the given `color` to the existing plot.
"""
function polyplot!(p, color)
    n = size(p, 2)
    if p[:,n] == p[:,1]
        plot!(p[1,:], p[2,:], linecolor=color)
    else
        plot!([p[1,:] p[1,1]], [p[2,:] p[2,1]], linecolor=color)
    end
end
