using LaTeXStrings
using PyPlot
using Printf

function dsdemo1()
    println("NTF Synthesis-- 5th-order modulator")

    fig = figure(1)
    fig.clf()
    fig.canvas.set_window_title("NTF Poles and Zeros")

    order = 5
    OSR = 32
    opt = 0
    H = synthesizeNTF(order, OSR, opt)

    plotPZ(H)

    print("paused")
    readline()

    fig = figure(2)
    fig.clf()
    fig.canvas.set_window_title("NTF Magnitude Response")
    fig.subplots_adjust(hspace=0.5)

    subplot(211)

    f = [range(0, 0.75/OSR, length=100); range(0.75/OSR, 0.5, length=100)]
    z = exp.(2im * pi * f)
    magH = dbv(evalTF(H, z))

    plot(f, magH, linewidth=1)
    figureMagic([0, 0.5], 0.05, 2, [-100, 10], 10, 2)
    xlabel(string("Normalize frequency ", L"(1 \rightarrow f_s)"))
    ylabel("dB")

    subplot(212)

    fstart = 0.01
    f = range(fstart, 1.2, length=200) / (2 * OSR)
    z = exp.(2im * pi * f)
    magH = dbv(evalTF(H, z))
    sigma_H = dbv(rmsGain(H, 0, 0.5/OSR))

    semilogx(f*2*OSR, magH, linewidth=1)
    semilogx([fstart, 1], sigma_H*[1, 1], linewidth=1, marker="o", markerfacecolor="none", markeredgewidth=0.5)
    text(0.15, sigma_H+2, @sprintf("rms gain = %5.0fdB", sigma_H))
    axis([fstart, 1.2, -100, -30])
    grid()
    grid(which="minor", linestyle="dotted")
    xlabel(string("Normalize frequency ", L"(1 \rightarrow f_B)"))
    ylabel("dB")

    println("paused")
    readline()


    println("Optimized zeros")

    fig = figure(1)
    fig.clf()
    fig.canvas.set_window_title("NTF Poles and Optimized Zeros")

    opt = 1
    H = synthesizeNTF(order, OSR, opt)

    plotPZ(H)

    fig.canvas.manager.window.raise_()

    print("paused")
    readline()

    fig = figure(2)
    axs = fig.get_axes()

    f = [range(0, 0.75/OSR, length=100); range(0.75/OSR, 0.5, length=100)]
    z = exp.(2im * pi * f)
    magH = dbv(evalTF(H, z))

    axs[1].plot(f, magH, linestyle="--", linewidth=1)

    f = range(fstart, 1.2, length=200) / (2 * OSR)
    z = exp.(2im * pi * f)
    magH = dbv(evalTF(H, z))

    axs[2].semilogx(f*2*OSR, magH, linestyle="--", linewidth=1)

    sigma_H = dbv(rmsGain(H, 0, 0.5/OSR))

    axs[2].plot([fstart, 1], sigma_H*[1, 1], linestyle="--", linewidth=1, marker="o", markerfacecolor="none", markeredgewidth=0.5)
    axs[2].text(0.15, sigma_H+2, @sprintf("rms gain = %5.0fdB", sigma_H))

    fig.canvas.manager.window.raise_()

    println("paused")
    readline()


    println("NTF Synthesis-- Bandpass Modulator")

    fig = figure(1)
    fig.clf()
    fig.canvas.set_window_title("Bandpass NTF Poles and Zeros")

    order = 8
    OSR = 64
    opt = 2
    f0 = 0.125
    H = synthesizeNTF(order, OSR, opt, 1.5, f0)

    plotPZ(H)

    fig.canvas.manager.window.raise_()

    print("paused")
    readline()

    fig = figure(2)
    fig.clf()
    fig.canvas.set_window_title("Bandpass NTF/STF Magnitude Response")

    subplot(211)

    f = [range(0, f0-1/(2*OSR), length=50); range(f0-1/(2*OSR), f0+1/(2*OSR), length=100); range(f0+1/(2*OSR), 0.5, length=50)]
    z = exp.(2im * pi * f)
    magH = dbv(evalTF(H, z))

    G = zpk(zeros(order÷2), H.matrix[1].p, 1, 1)
    G *= 1 / abs(evalTF(G, exp.(2im* pi * f0)))
    magG = dbv(evalTF(G, z))

    plot(f, magH, linewidth=1)
    plot(f, magG, color="r", linewidth=1)
    axis([0, 0.5, -100, 10])
    grid()
    xlabel(string("Normalize frequency ", L"(1 \rightarrow f_s)"))
    ylabel("dB")

    subplot(212)

    f = range(f0-0.3/OSR, f0+0.3/OSR, length=100)
    z = exp.(2im * pi * f)
    magH = dbv(evalTF(H, z))
    sigma_H = dbv(rmsGain(H, f0-0.25/OSR, f0+0.25/OSR))

    plot(2*OSR*(f .- f0), magH, linewidth=1)
    plot([-0.5, 0.5], sigma_H*[1, 1], marker="o", markerfacecolor="none", markeredgewidth=0.5)
    text(-0.2, sigma_H+2, @sprintf("rms gain = %5.0fdB", sigma_H))
    axis([-0.6, 0.6, -100, -60])
    grid()
    xlabel("Normalized frequency offset")
    ylabel("dB")

    fig.canvas.manager.window.raise_()

    return nothing
end
