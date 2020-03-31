using Plots
#using PyPlot
using ControlSystems

pyplot()



n = 3;
OSR = 64;
opt = 2;
Hinf = 1.7;
f0 = 1/16;
t = [0.5 1];
form = 'FB';
dbg = 1;

param, H, L0, ABCD = designLCBP(n,OSR,opt,Hinf,f0,t,form,x0,dbg)

#display(p8)
#display(v2)

