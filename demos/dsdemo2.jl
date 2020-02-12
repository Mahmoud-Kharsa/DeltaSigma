using ControlSystems
using DeltaSigma
using MAT
using Plots

OSR = 32
H = synthesizeNTF(5,OSR,1)
N = 8192
fB = ceil(N/(2*OSR))
ftest=floor(2/3*fB)
u = 0.5*sin.(2*pi*ftest/N*(0:N-1)')
v, = simulateDSM(u, H)

t = 0:100
plot(t, u[t .+ 1], line=:red, linetype=:steppre, legend=false)
plot!(t, v[t .+ 1], line=:green, linetype=:steppre, legend=false)
xlabel!("Sample Number")
ylabel!("u, v")