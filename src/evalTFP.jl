J = 2.0
b = 0.04
K = 1.0
R = 0.08
L = 1e-4

# Create the model transfer function
s = tf("s")
Hs = (K*s)/(s*((J*s + b)*(L*s + R) + K^2))


J2 = 2.0
b2 = 0.7
K2 = 1.0
R2 = 0.08
L2 = 1e-4

# Create the model transfer function
s2 = tf("s")
Hz = (K*s+2)/(s2*((J*s2 + b)*(L*s2 + R) + K^2))

f=[2;3;5]
function evalTFP(Hs,Hz,f)

    szeros, spoles, sk = zpkdata(Hs);
    szeros2 = szeros[1];
    spoles2 = spoles[1];
    zzeros, zpoles, zk = zpkdata(Hz);
    zzeros2 = zzeros[1];
    zpoles2 = zpoles[1];
    #return zpoles;
	slim = min(1e-3,max(1e-5,eps((1/(1+length(spoles))))));
    zlim = min(1e-3,max(1e-5,eps((1/(1+length(zzeros))))));

    H = zeros(size(f));
    w = 2*pi*f;	
    s = im*w;	
    z=exp.(s);
    #return zpoles; 
	for i=1:length(f)
        wi = w[i];
	    si = s[i];	
	    zi = z[i];
        if isempty(spoles2)
            cancel = false;
        else
            cancel = abs.(si-spoles2[i])<slim;
        end
        if !cancel
            #use evalRPoly directly with z/p/k data, much faster. (David Alldred, Feb 4 2013)
			#Only problem left is this line
            H[i] = sk .* evalRPoly(szeros2,si,1) ./ evalRPoly(spoles2,si,1) .* zk * evalRPoly(zzeros2,zi,1) ./ evalRPoly(zpoles2,zi,1);
        else
		    cancelz = abs.(zi-zzeros)<zlim;
		    if sum(cancelz) > sum(cancel)
		        H[i] = 0;
		    elseif sum(cancelz) < sum(cancel)
	    	    H[i] = Inf;
		    else
			    #Another problem of combination
			    H[i]=evalRPoly(szeros2,si,3) * zi^sum(cancel) * evalRPoly(zzeros2[!cancelz],zi,3) / (evalRPoly(spoles[!cancel],si)*evalRPoly(zpoles,zi,3))
		        #H[i] =  3;
		    end
        end

		    
    end
	return H
	end