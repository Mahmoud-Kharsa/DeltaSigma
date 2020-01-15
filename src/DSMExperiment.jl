using Plots
using PyPlot
using ControlSystems
using LinearAlgebra
#include("ds_quantize.jl")

function DSMExperiment(u,arg2, nlev, x0)
    nu = size(u,1)
    nq = length(nlev)
	ABCD=arg2
	order = size(ABCD,1)-nq
    #if(check_str(arg2))
    #    if (size(arg2,2) > 2) || (size(arg2,2)==nu+size(arg2,1))
	#        println("The ABCD argument does not have proper dimensions.");
	#	end
    #end
    #if isnan(x0)
    #x0 = zeros(order,1);
    #end
    #currently test with form=1, ABCD only.
    #p = zeros(order)
	#size(arg2,2) > 2 & size(arg2,2)==nu+size(arg2,1)
    #ABCD dimesions OK, example: size(arg2,2)=5, nu=2, size(arg2,1)=3
	#ABCD=arg2
	#order = size(ABCD,1)-nq;
	#if form_2 == 1
	A = ABCD[1:order, 1:order]
    B = ABCD[1:order, order+1:order+nu+nq]
    C = ABCD[order+1:order+nq, 1:order]
    D1= ABCD[order+1:order+nq, order+1:order+nu]
	#output:A=[3 1;4 3] B=[1 3;5 1] C=[1 2] D=[6]
	
	println("Matrix A: ", A)
	println("Matrix B: ", B)
	println("Matrix C: ", C)
	println("Matrix D: ", D1)
	#end
	global i
	N = length(u);
	v = zeros(nq,N);
    y = zeros(nq,N);
	#println(v)
	#println(N)
	#println(nq)
	#y = [[i] for i in 1.0:float(N)]
	#println(y)
	#y = C.*x0 + D1.*u
	#println(y[1:end])
	#y3=y[1:end]
	#println(y3[1])
	#for i=1:N
	#v[1]=ds_quantize(y[1],nq)
	#v[2]=ds_quantize(y[2],nq)
   # println(v[1])
	#println(v[2])
	#end
    for i=1:N
        y = C.*x0 + D1.*u
		#y = C.*x0 + D1.*u
		#print(y)

		#println(y[1])
        #y5, y6=ds_quantize(y[i],nlev);
		v[i]=ds_quantize(y[i],nq)
        #x0 = A .* x0 + B .* [u[:,i];v[:,i]];
		println("Quantizer Result")
		println(v)
		#i=i+1
	end
	#println(v[2])
end


DSMExperiment([2 7],[3 1 1 3 5;4 3 5 1 2;1 2 6 6 7],2,[7 12])
#DSMExperiment([2 5 6 4],[3 1 1;4 3 5;1 2 6],2,4)