using DeltaSigma
using Plots
using Printf

println("Invariant Set for MOD2 (find2dPIS)")

A = [1 0; 1 1]
B = [1 -1; 1 -2]
C = [0 1]
D = [0 0]
ABCD = [A B; C D]
order = 2

u = 1 / pi
s = find2dPIS(u, ABCD, 1)

N = 10000
skip = 100
_, x, = simulateDSM(u * ones(N+skip), ABCD, 2)
x = x[:,1+skip:N+skip]
nv = size(s, 2)
splus, eplus, sminus, eminus = dssplit2d(u, ABCD, s)
Buv = B * [u, 1]
s1 = A*splus + Buv*ones(1,size(splus, 2))
Buv = B * [u, -1]
s2 = A*sminus + Buv*ones(1,size(sminus, 2))
ns = [s1 s2]
out = outconvex2d(ns, s)

plot(legend=false, size=(600,500))

plot!(x[1,:], x[2,:], line=(:scatter), marker=(2, :black))

outi = (out .!= 0)[:]
plot!(ns[1,outi], ns[2,outi], line=(:scatter), marker=(:square, :white, stroke(:red)))

polyplot!(s, :blue)
polyplot!(s1, :magenta)
polyplot!(s2, :cyan)

str = @sprintf("Final Object: %d image vertices outside", sum(outi))
plot!(title=str)
plot!(xaxis=("x1", (-2, 2.5), -2:0.5:2.5), yaxis=("x2", (-3, 5), -3:1:5))
display(plot!())

@printf("%d points from the %d simulated states are outside.\n", sum(outconvex2d(x, s)), N)
@printf("%d image points are outside.\n", sum(out))
@printf("The returned polygon has %d vertices.\n", size(s, 2))
