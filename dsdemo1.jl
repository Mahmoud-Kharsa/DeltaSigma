using Plots
using PyPlot
using ControlSystems

pyplot()
H = tf([1], [1,2,1])
p = pzmap(H);
plot!(p, exp.(im*2*pi*(0:0.01:1)))