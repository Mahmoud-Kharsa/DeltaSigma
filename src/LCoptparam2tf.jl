function LCoptparam2tf(x,form)

    n = 3
    gw = ones(1,n-1)
    if form == "FB"
        gu = [1 zeros(1,n)];
        gv = x;
        gx = [zeros(1,n-1) 1];
        rx = zeros(1,n);
    elseif form ==  "FF"
        gu = [1 zeros(1,n-1) 1];
        gv = [1 zeros(1,n)];
        gx = x;
        rx = zeros(1,n);
    elseif form == "FFR"
        gu = [1 zeros(1,n-1) 1];
        gv = [x(1) zeros(1,n)];
        gx = [zeros(1,n-1) 1];
        rx = [0 x(2:n)];
    elseif form == "GEN"
        gv = [x(1:n)];
        rx = [0 x(n+1:2*n-1)];
    else
        @warn("The form listed is not supported.\n");
    end

#H, L0, ABCD = LCparam2tf(gu,gv,gx,rx,l,c);
    return H,Hpole,Hzero,gu,gw,gv,gx,rx,L0k
	end

