using ControlSystems
using LinearAlgebra
using Polynomials
function evalRPoly(roots,x,k)

    y = k;
    deleteat!(roots, findall(x->x==Inf, roots))
    for i in roots
        y = y*(x-i);
    end
	
	return y
	end
	
#roots2=[1, 2, 3];
#x=4;
#k=1;
#evalRPoly(roots2,x,k)