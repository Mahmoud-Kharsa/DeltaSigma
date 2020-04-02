using Printf

function dsdemo3()
    println("Modulator realization and scaling")
    println()

    fig = figure(1)
    fig.clf()
    fig.canvas.set_window_title("NTF")

    order = 5
    OSR = 42
    opt = 1
    H = synthesizeNTF(order, OSR, opt)
    a, g, b, c = realizeNTF(H)
    b = [b[1]; zeros(length(b)-1)]  # Use a single feed-in for the input

    plotPZ(H)

    println("Unscaled modulator")
    print("   DAC feedback coefficients = ")
    for i = 1:order
        @printf(" %.6f", a[i])
    end
    print("\n   resonator feedback coefficients = ")
    for i = 1:order÷2
        @printf(" %.6f", g[i])
    end
    println()
    println()

    println("Calculate the state maxima.")
    ABCD = stuffABCD(a, g, b, c)
    u = range(0, 0.6, length=30)
    N = 10000
    T = ones(N)
    maxima = zeros(order, length(u))
    for i = 1:length(u)
        ui = u[i]
        v, xn, xmax, = simulateDSM(ui * T, ABCD)
        maxima[:,i] = xmax[:]
        if any(xmax .> 1e2)
            umax = ui
            u = u[1:i]
            maxima = maxima[:,1:i]
            break
        end
    end

    fig = figure(2)
    fig.clf()
    fig.canvas.set_window_title("Simulated State Maxima")

    for i = 1:order
        semilogy(u, maxima[i,:], linestyle="--", linewidth=1, marker="o", markerfacecolor="none", markeredgewidth=0.5)
    end
    grid()
    grid(which="minor", linestyle="dotted")
    axis([0, 0.6, 1e-4, 10])
    xlabel("DC input")

    println("paused")
    readline()

    println("Calculate the scaled coefficients.")

    scaleABCD(ABCD, 2, 0, 1, NaN, NaN, 10000)
    ABCDs, umax, = scaleABCD(ABCD, 2, 0, 1, NaN, NaN, 10000)
    as, gs, bs, cs = mapABCD(ABCDs)

    @printf("Scaled modulator, umax=%.2f\n", umax)
    print("   DAC feedback coefficients = ")
    for i = 1:order
        @printf(" %.6f", as[i])
    end
    print("\n   resonator feedback coefficients = ")
    for i = 1:order÷2
        @printf(" %.6f", gs[i])
    end
    print("\n   interstage coefficients = ")
    for i = 1:order
        @printf(" %.6f", cs[i])
    end
    print("\n   feed-in coefficients = ")
    for i = 1:order
        @printf(" %.6f", bs[i])
    end
    println()
    println()

    println("Calculate the state maxima.")
    u = range(0, umax, length=30)
    N = 10000
    T = ones(N)
    maxima = zeros(order, length(u))
    for i = 1:length(u)
        ui = u[i]
        v, xn, xmax, = simulateDSM(ui * T, ABCDs)
        maxima[:,i] = xmax[:]
        if any(xmax .> 1e2)
            umax = ui
            u = u[1:i]
            maxima = maxima[:,1:i]
            break
        end
    end

    fig = figure(2)
    fig.clf()

    for i = 1:order
        semilogy(u, maxima[i,:], linestyle="--", linewidth=1, marker="o", markerfacecolor="none", markeredgewidth=0.5)
    end
    grid()
    xlabel("DC input")
    grid(which="minor", linestyle="dotted")
    axis([0, 0.6, 4e-2, 4])

    fig.canvas.manager.window.raise_()

    return nothing
end
