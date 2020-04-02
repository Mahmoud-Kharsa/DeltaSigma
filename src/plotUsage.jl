using PyPlot

"""
    plotUsage(sv, colors=["b","g","r","y"])

Plot the element usage for a multi-element DAC.
The colors are for `sv` = 1,-1,i,-i.
"""
function plotUsage(sv, colors=["b","g","r","y"])
    T = size(sv, 2)
    M = size(sv, 1)

    # Plot the grid
    x = repeat(0:T, inner=2)
    T2 = ceil(Int, (T+1)/2)
    y = [zeros(1, T2); M*ones(2, T2); zeros(1,T2)]
    y = y[1:2*(T+1)]
    plot(x, y, color="k")
    M2 = ceil(Int, (M+1)/2)
    x = [zeros(1, M2); T*ones(2, M2); zeros(1, M2)]
    x = x[1:2*(M+1)]
    y = repeat(0:M, inner=2)
    plot(x, y, color="k")
    axis("image")

    for t = 1:T
        for i = 1:M
            if sv[i,t] == 1
                fill([t-1, t-1, t, t], [i-1, i, i, i-1], color=colors[1])
            elseif sv[i,t] == -1
                fill([t-1, t-1, t, t], [i-1, i, i, i-1], color=colors[2])
            elseif sv[i,t] == 1im
                fill([t-1, t-1, t, t], [i-1, i, i, i-1], color=colors[3])
            elseif sv[i,t] == -1im
                fill([t-1, t-1, t, t], [i-1, i, i, i-1], color=colors[4])
            end
        end
    end
end
