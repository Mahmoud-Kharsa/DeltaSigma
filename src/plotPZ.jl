using ControlSystems
using PyPlot

"""
    plotPZ(H, color="b", markersize=5)

Plot the poles and zeros of a transfer function.
"""
function plotPZ(H, color="b", markersize=5)
    z, p, = zpkdata(H)

    # Plot x and o for poles and zeros, respectively
    plot(real(p[1]), imag(p[1]), color=color, linestyle="None", marker="x", markersize=markersize, markeredgewidth=0.5)
    plot(real(z[1]), imag(z[1]), color=color, linestyle="None", marker="o", markerfacecolor="none", markersize=markersize, markeredgewidth=0.5)

    # Draw unit circle, real axis and imag axis
    circle = exp.(im*2*pi*(0:0.01:1))
    plot(real(circle), imag(circle), linewidth=0.5, color=:darkorange)
    axis("equal")
    axis([-1, 1, -1, 1])
    plot([0, 0], [-2, 2], linestyle="dotted", linewidth=1, color=:black)
    plot([-2, 2], [0, 0], linestyle="dotted", linewidth=1, color=:black)
end
