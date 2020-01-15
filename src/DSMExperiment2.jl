using Plots
using PyPlot
using ControlSystems
using LinearAlgebra
using Printf
#include ("ds_quantize.jl")

#y=[-1 3 3]
#n=5
#global i
#for i=1:length(y)
    #y(:,i) = C*x0 + D1*u(:,i)
    #v(:,i) = ds_quantize(y(:,i),3)
	#show(v[1]); println()
#end
#function DSMExperiment2(u,arg2, nlev, x0)
#ABCD=arg2

struct Ntf
    poles
	zeros
end
ntf= Ntf([1],[1,2,1])

E = ss(zpk(ntf.poles,ntf.zeros,-1));	
order = length(ntf.zeros)
A=E.A
B2=E.B
C=E.C
D2=E.D
println("Original Matrix A", A)
println("Original Matrix B", B2)
println("Original Matrix C", C )
println("Original Matrix D", D2)
# Transform the realization so that C = [1 0 0 ...]
orderTemp=Matrix{Float64}(I, order, order)
Ctemp=[C' orderTemp]
F=qr(Ctemp)
Sinv = (F.Q)/norm(C); 
S = inv(Sinv);
C = C.*Sinv;
if C[1]<0
    S = -S;
	Sinv = -Sinv;
end
A = S.*Sinv.*(A); 
B2 = S.*B2; 
C = [1 zeros(1,order-1)]; 
D2 = 0;

B1 = -B2;
D1 = [1];
B = [B1 B2];
x0=[7;12;3]
u=[2;7]
println("Modified Matrix A", A)
println("Modified Matrix B", B)
println("Modified Matrix C", C)
println("Modified Matrix D", D1)
   # for i=1:length(u)
    #    y = C.*x0 + D1.*u
	#	#y = C.*x0
	#	print(y)

	#	println(y[1])
        #y5, y6=ds_quantize(y[i],nlev);
	#	v[i]=ds_quantize(y[i],nq)
        #x0 = A .* x0 + B .* [u[:,i];v[:,i]];
	#	println(v)
		#i=i+1
	#end 
  
#DSMExperiment2([2 7],[3 1 1 3 5;4 3 5 1 2;1 2 6 6 7],2,[7,12])