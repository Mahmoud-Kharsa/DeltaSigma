using Plots
using PyPlot
using ControlSystems
using LinearAlgebra
include("ds_optzeros.jl")

function designLCBP(n,OSR,opt,Hinf,f0,t,form,x0,dbg)
    optZeros=ds_optzeros(n,opt)
	if opt == 1
        w = 2*pi.*(f0 + 0.25/OSR.*optZeros');
    else
        w = 2*pi.*f0*ones(1,n);
    end
    l = (1/w);	
	c = (1/w);
	
	if x0 == NaN
        x = [2:n+1;]
    elseif length(x0) != n
        println("x0 has the wrong parameters.");
    else
        x = x0;
    end
	
	#if constr
	if dbg
		println("Iterations for initial pole placement.");
	end
	#options[2] = 0.001;	
	#options[3] = 1;	
	#options[14] = 1000;	
	#x = constr('LCObj1',x,options,[],[],[],param,0.97,dbg);
	H,Hpole,Hzero,gu,gw,gv,gx,rx,L0k= LCoptparam2tf(x,param);
	rmax = maximum(abs.(Hpole));
	if rmax>0.97
		println("Unable to find parameter values which stabilize the system.");
	end
	#if dbg>1
	 #   options[1] = 1; 	
	#end
	#if dbg == 1
	 #   println("\nParameter Optimization:\n");
	#end
	#options[2]  = 1e-4;	
	#options[3]  = 0.5;	
	#options[4]  = 0.01;	

	#options[14] = 1000;	
	#options[16] = 1e-4;	
	#options[17] = 1e-2;	
	#options[18] = .01;	 
	H,Hpole,Hzero,gu,gw,gv,gx,rx,L0k = LCoptparam2tf(x,param);
    H2=evalTFP(L0k,H,f0);

	gain_f0 = abs.(H2); 
	gu = gu/gain_f0; 
	L0k = L0k/gain_f0;
	ABCD[:,2*n+1] = ABCD[:,2*n+1]/gain_f0;
	if form=="FF" || form =="FFR"
		gu[n+1] =1;
	end

    #if dbg ==1
	 #   LCplotTF(H,L0,param);
    #end
	return param, H, L0, ABCD
	end