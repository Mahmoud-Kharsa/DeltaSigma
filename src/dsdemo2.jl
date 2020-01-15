using Plots
#using PyPlot
using ControlSystems

pyplot()



#include("synthesizeNTF.jl")
OSR = 32;
H = synthesizeNTF(5,OSR,1);
#N = 8192;
N = 100;
fB = ceil(N/(2*OSR)); 
ftest=floor(2/3*fB);
u = 0.5*sin.(2*pi*ftest/N*(0:N-1));	
#v = simulateDSM([2 7],[3 1 1 3 5;4 3 5 1 2;1 2 6 6 7],2,[7,12]); 
#v = 3;
u=transpose(u)
#v=transpose(v)
p8=plot(0:N-1, u, line = :red, linetype=:steppre, title="dsdemo2_1st_example", xlabel="Sample Number", ylabel="value of u and v")
#v2=plot(0:N-1, v, line = :green, linetype=:steppre, title="dsdemo2_1st_example", xlabel="Sample Number", ylabel="value of u and v")

display(p8)
#display(v2)

