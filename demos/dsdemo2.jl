using ControlSystems
using DeltaSigma
using Plots

OSR = 32
H = synthesizeNTF(5,OSR,1)
N = 8192
fB = ceil(N/(2*OSR))
ftest = floor(2/3*fB)
u = 0.5*sin.(2*pi*ftest/N*(0:N-1)')
v, = simulateDSM(u, H)

t = 0:100
p1 = plot(t, u[t .+ 1], line=:red, linetype=:steppre, legend=false)
plot!(p1, t, v[t .+ 1], line=:green, linetype=:steppre, legend=false)
xlabel!("Sample Number")
ylabel!("u, v")
display(p1)

snr, amp = simulateSNR(H, OSR)
p4 = plot(amp[:], snr[:], seriestype=:scatter, legend=false)
plot!(p4, xlims=(-100,0), xticks=-100:10:0)
plot!(p4, ylims=(0, 100), yticks=0:10:100)
display(p4)
