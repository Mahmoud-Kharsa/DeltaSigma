using Polynomials
function LCparam2tf(gu,gv,gx,rx, l, cond)

    f0 = mean(1 ./ sqrt.(l.*cond)) ./ (2*pi);
    n = length(l);

    if isempty(k) || k==NaN
        k = 1;
    end


    t = [0 1];


    if length(rl)==1	
	    rl = [ones(1,n)];

    else
        rl = [zeros(1,n)];
    end


	# Form the state-space description of the modulator
	# The states are ordered v1, i1, v2, i2 ...
	n2 = 2*n;
	ABc = zeros(n2,n2+2);
	CDc = zeros(1,n2+2);
	for i = 1:n
		y1=[1 zeros(1,n2+1)]
		i2 = 2*i-1;
		# y represents the voltage at the top of the resistor stacked on the tank.
		#need to adjust calculation later, no problem on functionality. 
		if i == 1
			y = [ 1 zeros(1,n2+1) ];
		else
			ABc[i2,:] = gw[i-1].*y1;
			y = rx[i].*gw[i-1].*y1; 
			y[i2] = 1;
		end
		ABc[i2,n2+1:n2+2] = [gu[i] -1*gv[i]]; 
		ABc[i2:i2+1,i2:i2+1] = [-gc[i] -1; 1 -rl[i]];
		ABc[i2,:] = ABc[i2,:]./cond[i]; 
		ABc[i2+1,:] = ABc[i2+1,:]./l[i]; 
		CDc2=zeros(1,n2+2);
		CDc = CDc2 + gx[i].*y1;
	end
	CDc = CDc + [zeros(1,n2) gu(n+1) -1*gv(n+1)];


	Ac,Bc,Cc,Dc = partitionABCD([ABc;CDc], 2);
	#sys_c = ss( Ac, Bc[:,2], Cc, Dc[2] );
	#sys = mapCtoD(sys_c,t,f0);
	# Augment the input matrix. 
	#A=sys.a;
	#B=[padb(Bc(:,1),size(sys.b,1)) sys.b];
	#C=sys.c;
	#D=[Dc(1) sys.d];
	A=Ac
	B=Bc[:,2]
	C=Cc
	D=Dc[2]

	#Compute L0; use the LC parameters to compute the poles and thereby
	# ensure exact agreement with the zeros in H.
	s_poles = zeros(1,2*n);

	for i=1:n
		s_poles[2*i-1:2*i] = roots(Poly([1+rl[i].*gc[i],( rl[i].*cond[i] + gc[i].*l[i] ),l[i].*cond[i]]));
	end
	LF0 = ss(Ac,Bc[:,1],Cc,Dc[1]);
	L0z, L0p, L0k = zpkdata( LF0 );
	L0p =  s_poles;

	# Compute H. Use the poles in L0 to compute the zeros in H.
	ABCD =[A B; C D];
	#if k==0		
	#	w0 = mean(1./sqrt(l.*cond));
	#	H = calculateTF(ABCD);
	#	Hzero = exp(s_poles);
	#	stf0 = abs(evalTFP(L0,H,w0/(2*pi)));
	#	u = 0.1/stf0*sin(w0*[0:10000]);
	#	[tmp1,tmp2,tmp3,y] = simulateDSM(u,[A B; C D]);
	#	k = mean(abs.(y))/mean(y.^2);
	#end
	sys1=ss(A,B,C,D,k)
	H = tf(sys1);
	HZero = exp.(s_poles);

	# Correct L0k to include the quantizer gain
	L0k = L0k.*k;
	
	return H, ABCD, L0k
	end