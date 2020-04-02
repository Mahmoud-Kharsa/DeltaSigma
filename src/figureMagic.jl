using PyPlot

function figureMagic(xRange, dx, xLab, yRange, dy, yLab)
    axis([xRange; yRange])
    grid(b=true)

    xtix = collect(xRange[1]:dx:xRange[2])
    xticks(xtix, axisLabels(xtix, xLab))

    ytix = collect(yRange[1]:dy:yRange[2])
    yticks(ytix, axisLabels(ytix, yLab))
end
