using FFTW
using LaTeXStrings
using Printf
using PyPlot

function dsdemo2()
    println("Discrete-Time Simulation")

    fig = figure(1)
    fig.clf()
    fig.canvas.set_window_title("Modulator Input & Output")

    OSR = 32
    H = synthesizeNTF(5, OSR, 1)
    N = 8192
    fB = ceil(Int, N/(2*OSR));
    ftest = floor(Int, 2/3*fB);
    u = 0.5 * sin.(2*pi*ftest/N * (0:N-1))
    v, = simulateDSM(u, H)
    t = 0:100

    step(t, u[t.+1], color="r", linewidth=1)
    step(t, v[t.+1], color="g", linewidth=1)
    axis([0, 100, -1.2, 1.2])
    xlabel("Sample Number")
    ylabel("u, v")

    print("paused")
    readline()

    fig = figure(2)
    fig.clf()
    fig.canvas.set_window_title("Output Spectrum")

    f = range(0, 0.5; length=N÷2+1)
    spec = fft(v .* ds_hann(N) / (N/4))

    plot(f, dbv(spec[1:N÷2+1]), color="b", linewidth=1)
    figureMagic([0, 0.5], 0.05, 2, [-120, 0], 20, 1)
    xlabel("Normalized Frequency")
    ylabel("dBFS")

    print("paused")
    readline()

    snr = calculateSNR(spec[3:fB+1], ftest-2)

    text(0.05, -10, @sprintf("SNR = %4.1fdB @ OSR = %d", snr, OSR))

    print("paused")
    readline()

    NBW = 1.5 / N
    Sqq = 4 * evalTF(H, exp.(2*im*pi*f)).^2 / 3

    plot(f, dbp(Sqq*NBW), color="m", linewidth=2)
    text(0.5, -90, string(@sprintf("NBW = %4.1E × ", NBW), L"f_s"), ha="right")
    legend(["Simulation", "Expected PSD"])

    print("paused")
    readline()

    fig = figure(3)
    fig.clf()
    fig.canvas.set_window_title("SQNR")

    snr_pred, amp_pred = predictSNR(H, OSR)
    snr, amp = simulateSNR(H, OSR)
    pk_snr, pk_amp = peakSNR(snr, amp)

    plot(amp_pred, snr_pred, linewidth=1)
    plot(amp, snr, color="g", linestyle="None", marker="o", markerfacecolor="none", markeredgewidth=0.5)
    text(-25, 85, @sprintf("peak SNR = %4.1fdB\n@ OSR = %d\n", pk_snr, OSR), ha="right")
    figureMagic([-100, 0], 10, 1, [0, 100], 10, 1)
    xlabel("Input Level (dBFS)")
    ylabel("SQNR (dB)")

    println("paused")
    readline()


    println("Bandpass Modulator")

    fig = figure(1)
    fig.clf()

    f0 = 1/8
    OSR = 64
    H = synthesizeNTF(8, OSR, 1, 1.5, f0)
    fB = ceil(Int, N/(2*OSR))
    ftest = round(Int, f0*N + 1/3*fB)
    u = 0.5 * sin.(2*pi*ftest/N * (0:N-1))   # half-scale sine-wave input
    v, = simulateDSM(u, H)
    t = 0:100

    step(t, u[t.+1], color="r", linewidth=1)
    step(t, v[t.+1], color="g", linewidth=1)
    axis([0, 100, -1.2, 1.2])
    xlabel("Sample Number")
    ylabel("u, v")

    fig.canvas.manager.window.raise_()

    print("paused")
    readline()

    fig = figure(2)
    fig.clf()

    f = range(0, 0.5; length=N÷2+1)
    spec = fft(v .* ds_hann(N) / (N/4))
    f1 = round(Int, (f0-0.25/OSR)*N)
    f2 = round(Int, (f0+0.25/OSR)*N)
    snr = calculateSNR(spec[f1:f2], ftest-f1+1)

    plot(f, dbv(spec[1:N÷2+1]), color="b", linewidth=1)
    text(0.15, -10, @sprintf("SNR = %4.1fdB @ OSR = %d", snr, OSR))
    figureMagic([0, 0.5], 0.05, 2, [-140, 0], 20, 1)
    xlabel("Normalized Frequency")
    ylabel("dBFS")

    fig.canvas.manager.window.raise_()

    print("paused")
    readline()

    NBW = 1.5/N
    Sqq = 4 * evalTF(H, exp.(2*im*pi*f)).^2 / 3

    plot(f, dbp(Sqq*NBW), color="m", linewidth=2)
    text(0.5, -90, string(@sprintf("NBW = %4.1E × ", NBW), L"f_s"), ha="right")
    legend(["Simulation", "Expected PSD"])

    print("paused")
    readline()

    fig = figure(3)
    fig.clf()

    snr_pred, amp_pred = predictSNR(H, OSR, [-120:10:-20; -15; -10:0], f0)
    snr, amp = simulateSNR(H, OSR, [-120:10:-20; -15; -10:0], f0)
    pk_snr, pk_amp = peakSNR(snr, amp)

    plot(amp_pred, snr_pred, linewidth=1)
    plot(amp, snr, color="g", linestyle="None", marker="o", markerfacecolor="none", markeredgewidth=0.5)
    text(-20, 95, @sprintf("peak SNR = %4.1fdB\n@ OSR = %d\n", pk_snr, OSR), ha="right")
    figureMagic([-110, 0], 10, 1, [0, 110], 10, 1)
    xlabel("Input Level (dBFS)")
    ylabel("SQNR (dB)")

    fig.canvas.manager.window.raise_()

    println("paused")
    readline()

    fig = figure(3)
    fig.clf()
    fig = figure(2)
    fig.clf()
    fig = figure(1)
    fig.clear()
    fig.canvas.set_window_title("Input & Output")
    fig.subplots_adjust(hspace=0.5)

    colors = ["m", "b"]
    Hinf_list = [2, 8]
    for i = 1:2
        Hinf = Hinf_list[i]
        col = colors[i]
        println("15-step 7th-order Lowpass")
        println("Hinf = $Hinf")

        fig = figure(1)

        OSR = 8
        M = 16
        H = synthesizeNTF(7, OSR, 1, Hinf)
        N = 8192
        fB = ceil(Int, N/(2*OSR))
        ftest = floor(Int, 2/7*fB)
        u = 0.5*M*sin.(2*pi*ftest/N*(0:N-1))
        v, = simulateDSM(u, H, M+1)
        t = 0:100

        subplot(2, 1, i)

        step(t, u[t.+1], color="g", linewidth=1)
        step(t, v[t.+1], color=col, linewidth=1)
        figureMagic([0, 100], 10, 2, [-M, M], 2, 4)
        xlabel("Sample Number")
        ylabel("u, v")

        fig.canvas.manager.window.raise_()

        print("paused")
        readline()

        fig = figure(2)

        f = range(0, 0.5; length=N÷2+1)
        spec = fft(v .* ds_hann(N) / (M*N/4))
        snr = calculateSNR(spec[3:fB+1], ftest-2)
        NBW = 1.5/N
        Sqq = 4 * evalTF(H, exp.(2*im*pi*f)).^2 / (3*M^2)

        plot(f, dbv(spec[1:N÷2+1]), color=col)
        plot(f, dbp(Sqq*NBW), color="c", linewidth=2)
        text(0.1, 10*(i-3), @sprintf("SNR = %4.1fdB @ OSR = %d", snr, OSR), color=col)
        if i == 1
            text(0.5, -110, string(@sprintf("NBW = %4.1E × ", NBW), L"f_s"), ha="right")
        end
        figureMagic([0, 0.5], 0.05, 2, [-160, 0], 20, 1)
        xlabel("Normalized Frequency")
        ylabel("dBFS")

        fig.canvas.manager.window.raise_()

        println("paused")
        readline()

        fig = figure(3)

        snr, amp = simulateSNR(H, OSR, vcat(-120:10:-20, -15, -10:0), 0, M+1)
        pk_snr, pk_amp = peakSNR(snr, amp)

        plot(amp, snr, color=col, linestyle="None", marker="o", markerfacecolor="none", markeredgewidth=0.5)
        plot(amp, snr, color=col, linestyle="--", linewidth=1)
        text(-13, pk_snr, @sprintf("peak SNR = %4.1fdB\n@ OSR = %d\n", pk_snr, OSR), fontsize=8, color=col, ha="right", va="top")
        figureMagic([-120, 0], 10, 2, [0, 120], 10, 2)
        xlabel("Input Level (dBFS)")
        ylabel("SQNR (dB)")

        fig.canvas.manager.window.raise_()
    end

    return nothing
end
