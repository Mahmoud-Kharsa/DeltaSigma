using Plots
using ControlSystems

pyplot()


#Plot Poles & Zeros 1
order = 5
OSR = 32
opt = 0
H = synthesizeNTF(order,OSR,opt)
H1 = H

p1 = plot(overwrite_figure = true);
plot!(size = (650, 530), aspect_ratio = 1)
pzmap!(H1);
plot!(exp.(im*2*pi*(0:0.01:1)))

display(p1)
readline()


#Plot Magnitude 1
f = [collect(range(0,0.75/OSR,length=100)); collect(range(0.75/OSR,0.5,length=100))];
z = exp.(2im * pi * f);
magH = H(z,false)[:,1,1];
magH = 20*log10.(abs.(magH));

fstart = 0.01;
f2 = collect(range(fstart, 1.2, length=200))/(2 * OSR); 
z2 = exp.(2im * pi * f2);
magH2 = H(z2,false)[:,1,1];
magH2 = 20*log10.(abs.(magH2));

p4 = plot(f, magH, overwrite_figure = true);
plot!(xlims = (0,0.5), ylims = (-100,10))
plot!(size = (700, 530))

p5 = plot(f2*2*32, magH2, overwrite_figure = true);
plot!([0.01; 1], [-47.365252553595320; -47.365252553595320])
scatter!([0.01; 1], [-47.365252553595320; -47.365252553595320])
plot!(xlims = (0.01,1.2), ylims = (-100,-30), xaxis=:log10)
plot!(size = (700, 530))

p6 = plot(p4,p5,layout=(2,1),legend=false,overwrite_figure = true);

display(p6)
readline()


#Plot Poles & Zeros 2
opt = 1;
H = synthesizeNTF(order,OSR,opt);
H2 = H;

p2 = plot(overwrite_figure = true);
plot!(size = (650, 530), aspect_ratio = 1)
pzmap!(H2);
plot!(exp.(im*2*pi*(0:0.01:1)))


display(p2)
readline()

#Plot Poles & Zeros 3
order = 8;
OSR = 64;
opt = 2;
f0 = 0.125;
H = synthesizeNTF(order,OSR,opt, 1.5 ,f0);
H3 = H;

p3 = plot(overwrite_figure = true);
plot!(size = (650, 530), aspect_ratio = 1)
pzmap!(H3);
plot!(exp.(im*2*pi*(0:0.01:1)))

display(p1)
