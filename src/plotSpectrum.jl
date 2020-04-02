using PyPlot

"""
    plotSpectrum(X, fin, col, nbin=8, n=3)

Plot a smoothed spectrum
"""
function plotSpectrum(X, fin, col, nbin=8, n=3)
    f, p = logsmooth(X, fin, nbin, n)
    semilogx(f, p, color=col, linewidth=0.5)
end
