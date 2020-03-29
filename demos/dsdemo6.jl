using DeltaSigma
using DSP
using FFTW
using Plots

println("Half-band filter design")
println("ALPHA VERSION")

if true
    # Use one of the canned examples
    f1, f2 = exampleHBF(2)
else
    # designHBF
end

# interleave the even and odd decimated impulse responses
Nimp = 2^11
imp = simulateHBF([1; zeros(Nimp-1)], f1[1], f2[1])

mag = abs.(fft(imp))
mag = mag[1:end÷2+1]
plt = plot(legend=false, size=(550,450))
plot!(plt, range(0, 0.5, length=length(mag)), amp2db.(abs.(mag)))
plot!(plt, xlims=(0, 0.5), ylims=(-150, 3))
plot!(plt, xticks=(0:0.05:0.5, ["0", "", "0.1", "", "0.2", "", "0.3", "", "0.4", "", "0.5"]))
plot!(plt, yticks=(-150:10:3, ["-150", "", "", "", "", "-100", "", "", "", "", "-50", "", "", "", "", "0"]))
display(plt)
