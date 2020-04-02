using FFTW
using PyPlot

function dsdemo6()
    println("Half-band filter design")
    println("ALPHA VERSION")

    if false
        # Use one of the canned examples
        f1, f2 = exampleHBF(2)
    else
        fp = 0.9 * 0.25
        delta = undbv(-100)
        fig = figure(1)
        fig.clf()
        fig.canvas.set_window_title("designHBF Iterations")
        f1, f2, = designHBF(fp, delta, true)
    end

    # interleave the even and odd decimated impulse responses
    Nimp = 2^11
    imp = simulateHBF([1; zeros(Nimp-1)], f1[1], f2[1])

    mag = abs.(fft(imp))
    mag = mag[1:end÷2+1]

    fig = figure(2)
    fig.clf()
    fig.canvas.set_window_title("designHBF Result")

    plot(range(0, 0.5, length=length(mag)), dbv(mag))
    figureMagic([0, 0.5], 0.05, 2, [-150, 3], 10, 5)

    return nothing
end
