using PyPlot

"""
    polyplot(p, color)

Plot a polygon given by the point list `p` with the given `color`.
"""
function polyplot(p, color)
    n = size(p, 2)
    if p[:,n] == p[:,1]
        plot(p[1,:], p[2,:], color=color, linewidth=1)
    else
        plot([p[1,:] p[1,1]], [p[2,:] p[2,1]], color=color, linewidth=1)
    end
end
