using ControlSystems
using DeltaSigma
using DSP
using FFTW
using Plots
using Printf

println("Discrete-Time Simulation")

OSR = 32
H = synthesizeNTF(5, OSR, 1)
N = 8192
fB = ceil(Int, N/(2*OSR));
ftest = floor(Int, 2/3*fB);
u = 0.5 * sin.(2*pi*ftest/N*(0:N-1))
v, = simulateDSM(u, H)

t = 0:100
p1 = plot(legend=false, size=(590,470))
plot!(p1, t, u[t .+ 1], line=(:steppre, :red))
plot!(p1, t, v[t .+ 1], line=(:steppre, :green1))
plot!(p1, xaxis=("Sample Number", (0,100), 0:10:100), yaxis=("u, v", (-1.2, 1.2)))
display(p1)
print("paused")
readline()

f = range(0, 0.5; length=N÷2+1)
spec = fft(v .* hanning(N+1)[1:N]/(N/4))
p2 = plot(legend=false, size=(590,470))
plot!(p2, f, amp2db.(abs.(spec[1:N÷2+1])), label="Simulation", line=:blue)
plot!(p2, xaxis=("Normalized Frequency", (0,0.5)), yaxis=("dbFS", (-120,0)))
plot!(p2, xticks=(0:0.05:0.5, ["0", "", "0.1", "", "0.2", "", "0.3", "", "0.4", "", "0.5"]))
display(p2)
print("paused")
readline()

snr = calculateSNR(spec[3:fB+1], ftest-2)
plot!(p2, ann=(0.05, -10, text(@sprintf("SNR = %4.1fdB @ OSR = %d", snr, OSR), 10, :left)))
display(p2)
print("paused")
readline()

NBW = 1.5/N
Sqq = 4 * H(exp.(2*im*pi*f), false)[:,1,1].^2 / 3
plot!(p2, f, pow2db.(abs.(Sqq*NBW)), label="Expected PSD", line=(2, :magenta), legend=true)
plot!(p2, ann=(0.5, -90, text(@sprintf("NBW = %4.1E x fs", NBW), 10, :right)))
display(p2)
print("paused")
readline()

snr_pred, amp_pred = predictSNR(H, OSR)
snr, amp = simulateSNR(H, OSR)
pk_snr, pk_amp = peakSNR(snr, amp)
p3 = plot(legend=false, size=(590,470))
plot!(amp_pred, snr_pred)
plot!(p3, amp, snr, line=(:scatter), marker=(5, :white, stroke(:green)))
plot!(p3, xaxis=("Input Level (dBFS)", (-100,0), -100:10:0), yaxis=("SQNR (dB)", (0, 100), 0:10:100))
plot!(p3, ann=(-25, 85, text(@sprintf("peak SNR = %4.1fdB\n@ OSR = %d\n", pk_snr, OSR), 10, :right)))
display(p3)
print("paused")
readline()

println("Bandpass Modulator")

f0 = 1/8
OSR = 64
H = synthesizeNTF(8, OSR, 1, 1.5, f0)
fB = ceil(Int, N/(2*OSR))
ftest = round(Int, f0*N + 1/3*fB)
u = 0.5*sin.(2*pi*ftest/N*(0:N-1))   # half-scale sine-wave input
v, = simulateDSM(u, H)

t = 0:100
p1 = plot(legend=false, size=(590,470))
plot!(p1, t, u[t .+ 1], line=(:steppre, :red))
plot!(p1, t, v[t .+ 1], line=(:steppre, :green1))
plot!(p1, xaxis=("Sample Number", (0,100), 0:10:100), yaxis=("u, v", (-1.2, 1.2)))
display(p1)
print("paused")
readline()

f = range(0, 0.5; length=N÷2+1)
spec = fft(v .* hanning(N+1)[1:N]/(N/4))
f1 = round(Int, (f0-0.25/OSR)*N)
f2 = round(Int, (f0+0.25/OSR)*N)
snr = calculateSNR(spec[f1:f2], ftest-f1+1)
p2 = plot(legend=false, size=(590,470))
plot!(p2, f, amp2db.(abs.(spec[1:N÷2+1])), label="Simulation", line=:blue)
plot!(p2, xaxis=("Normalized Frequency", (0,0.5), 0:0.05:0.5), yaxis=("dBFS/NBW", (-140,0)))
plot!(p2, ann=(0.15, -10, text(@sprintf("SNR = %4.1fdB @ OSR = %d", snr, OSR), 10, :left)))
display(p2)
print("paused")
readline()

NBW = 1.5 / N
Sqq = 4 * H(exp.(2*im*pi*f), false)[:,1,1].^2 / 3
plot!(p2, f, pow2db.(abs.(Sqq*NBW)), label="Expected PSD", line=(2, :magenta), legend=true)
plot!(p2, ann=(0.5, -90, text(@sprintf("NBW = %4.1E x fs", NBW), 10, :right)))
display(p2)
print("paused")
readline()

snr_pred, amp_pred = predictSNR(H, OSR, vcat(-120:10:-20, -15, -10:0), f0)
snr, amp = simulateSNR(H, OSR, vcat(-120:10:-20, -15, -10:0), f0)
pk_snr, pk_amp = peakSNR(snr, amp)
p3 = plot(legend=false, size=(590,470))
plot!(p3, amp_pred, snr_pred)
plot!(p3, amp, snr, line=(:scatter), marker=(5, :white, stroke(:green)))
plot!(p3, xaxis=("Input Level (dBFS)", (-110,0), -110:10:0), yaxis=("SQNR (dB)", (0, 110), 0:10:110))
plot!(p3, ann=(-20, 95, text(@sprintf("peak SNR = %4.1fdB\n@ OSR = %d\n", pk_snr, OSR), 10, :right)))
display(p3)
print("paused")
readline()

colors = [:magenta, :blue]
Hinf_list = [2, 8]
p1 = plot(layout=(2,1), legend=false, size=(590,470))
p2 = plot(legend=false, size=(590,470))
p3 = plot(legend=false, size=(590,470))
for i = 1:2
    Hinf = Hinf_list[i]
    col = colors[i]
    println("15-step 7th-order Lowpass")
    println("Hinf = $Hinf")
    OSR = 8
    M = 16
    H = synthesizeNTF(7, OSR, 1, Hinf)
    N = 8192
    fB = ceil(Int, N/(2*OSR))
    ftest = floor(Int, 2/7*fB)
    u = 0.5*M*sin.(2*pi*ftest/N*(0:N-1))
    v, = simulateDSM(u, H, M+1)

    plot!(p1, t, u[t .+ 1], line=(:steppre, :green1), subplot=i)
    plot!(p1, t, v[t .+ 1], line=(:steppre, col), subplot=i)
    plot!(p1, xaxis=("Sample Number", (0,100), 0:10:100), yaxis=("u, v", (-M, M)))
    display(p1)
    print("paused")
    readline()

    f = range(0, 0.5; length=N÷2+1)
    spec = fft(v .* hanning(N+1)[1:N]/(M*N/4))
    snr = calculateSNR(spec[3:fB+1], ftest-2)
    NBW = 1.5/N
    Sqq = 4 * H(exp.(2*im*pi*f), false)[:,1,1].^2 / (3*M^2)
    plot!(p2, f, amp2db.(abs.(spec[1:N÷2+1])), line=col)
    plot!(p2, xaxis=("Normalized Frequency", (0,0.5)), yaxis=("dBFS", (-160,0)))
    plot!(p2, xticks=(0:0.05:0.5, ["0", "", "0.1", "", "0.2", "", "0.3", "", "0.4", "", "0.5"]))
    plot!(p2, ann=(0.1, 10*(i-3), text(@sprintf("SNR = %4.1fdB @ OSR = %d", snr, OSR), col, 10, :left)))
    plot!(p2, f, pow2db.(abs.(Sqq*NBW)), line=(2, :cyan))
    plot!(p2, ann=(0.5, -110, text(@sprintf("NBW = %4.1E x fs", NBW), 10, :right)))
    display(p2)
    print("paused")
    readline()

    snr, amp = simulateSNR(H, OSR, vcat(-120:10:-20, -15, -10:0), 0, M+1)
    pk_snr, pk_amp = peakSNR(snr, amp)
    plot!(p3, amp, snr, line=(:dash, col), marker=(5, :white, stroke(col)))
    plot!(p3, xaxis=("Input Level (dBFS)", (-120,0), -120:10:0), yaxis=("SNR (dB)", (0, 120), 0:10:120))
    plot!(p3, ann=(-13, pk_snr, text(@sprintf("peak SNR = %4.1fdB\n@ OSR = %d\n", pk_snr, OSR), col, 10, :right)))
    display(p3)
end

nothing
