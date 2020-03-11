using Pkg
Pkg.add("Roots")
using Roots
Pkg.add("Polynomials")
using Polynomials

function realizeNTF(ntf,form,stf)
    #function returns a transfer function
    # utilizes the return function of synthesizedNTF to feed into realize NTF
    # form=> only implement CRFB, CRFF, CIFB, CIFF
    # a => Feedback/feedforward coefficients from/to the quantizer (1 × n)
    # g => Resonator coefficients (1 × ⌊n/2⌋)
    # b => Feed-in coefficients from the modulator input to each integrator (1×(n+1))
    # c => Integrator inter-stage coefficients. (1 × n). In unscaled modulators, c is all ones

    # # [a,g,b,c] = realizeNTF(ntf,form='CRFB',stf=1)
    # # Convert a noise transfer function into coefficients for the desired structure.
    # # Supported structures are
    # #	CRFB	Cascade of resonators, feedback form.
    # # 	CRFF	Cascade of resonators, feedforward form.
    # #	CIFB	Cascade of integrators, feedback form.
    # #	CIFF	Cascade of integrators, feedforward form.
    # # 	CRFBD	CRFB with delaying quantizer.
    # # 	CRFFD	CRFF with delaying quantizer.
    # #   PFF     Parallel feed-forward.
    # #	Stratos A CIFF-like structure with non-delaying resonator feedbacks,
    # #               contributed in 2007 by Jeff Gealow
    # #   DSFB    Cascade of double-sampled integrators, feedback form.
    # # See the accompanying documentation for block diagrams of each structure
    # #
    # # The order of the NTF zeros must be (real, complex conj. pairs).
    # # The order of the zeros is used when mapping the NTF onto the chosen topology.
    # #
    # # stf is a zpk transfer function

    # # The basic idea is to equate the loop filter at a set of
    # # points in the z-plane to L1 = 1-1/ntf at those points.
    # stderr = 2
    # # Handle the input arguments

    parameters = ["ntf","form","stf"]
    #NaN stands for not a number
    defaults = [NaN, "CRFB", []]
    # for i=1:length(defaults)
    # #length of default should be 3
    #     parameter = char(parameters(i))
    #     if i>nargin | ( eval(['isnumeric(' parameter ') '])  &  ...
    #             eval(['any(isnan(' parameter ')) | isempty(' parameter ') ']) )
    #         eval([parameter '=defaults{i}'])
    #     # nargin returns the number of input arguments used to call the function
    #     # make sure it doesn't exceed 4
    #     # true aka 1 if any number element of a vector is non-zero
    #     # true aka 1 if not a number
    #     end
    # end
    # form = "CRFB"
    ntf_z = [1.00000000000000 + 0.00000000000000im;0.999188974764100 + 0.0402665209555813im;0.999188974764100 - 0.0402665209555813im;0.997703673259689 + 0.0677302027468142im;0.997703673259689 - 0.0677302027468142im]
    ntf_p = [0.807058252221672 - 0.119565829763585im 0.807058252221672 + 0.119565829763585im 0.898412456856058 - 0.219273016428406im 0.898412456856058 + 0.219273016428406im 0.778301510991778 + 0.00000000000000im]
    # ntf_p_roots = roots([0.807058252221672 - 0.119565829763585im, 0.807058252221672 + 0.119565829763585im, 0.898412456856058 - 0.219273016428406im, 0.898412456856058 + 0.219273016428406im, 0.778301510991778 + 0.00000000000000im])
     # Code common to all functions
    # ntf_p = ntf.p{1}
    # within the return object H of Synthesize NTF, P contains poles, and the at 1,1 of the H.P
    # contains real and imaginary values stored in an array
    # ntf_z = ntf.z{1}
    # applies the same for zeroes
    # println("nt_zf",ntf_z_roots)
    order = length(ntf_p)
    # the order is calculated by evaluating the length of the pole vector
    order2 = order ÷ 2
    # order 2 by rounding down after division of 2
    odd = order - 2*order2
    # process odd order
    # zeroes function create a M x N matrix of 0s, M,N
    a = zeros(Float64, 1, order)
    # a = zeros(1,order)
    g = zeros(Float64, 1, order2)
    # g = zeros(1,order2)
    b = zeros(Float64, 1, order+1)
    # b = zeros(1,order+1)
    c = ones(Float64, 1, order)
    # ones function create a M x N matrix of 1s, M,N
    # c = ones(1,order)
    # basically initializing matrices for the return function
    stf=[]
     # Choose a set of points in the z-plane at which to try to make L1 = 1-1/H
     # I don't know how to do this properly, but the code below seems to work
     # for order <= 11
    N = 200
    min_distance = 0.09
    # C is a matrix of N by 1 of 0s
    C = zeros(ComplexF64, N, 1)
    # C = zeroes(N,1)
    j = 1
    for i=1:N
        # do 200 iterations
        z = 1.1*exp(2i*pi*im/N)
        # println("z value: ",z)
        # ntf_z is the zeroes of real and imaginary extracted from the H function, what is z
        if all(ntf_z->abs(ntf_z-z)>min_distance,ntf_z)
            # jth element as z
            C[j] = z 
            # z calculates some value and it looks like it assigns the value
            # at a particular j 
            j = j+1
            # then it proceeds to increment
        end
    end
    # println("C Matrix",C)
    # for the rest of the j that do not meet the criteria assign null vector
    # for loop from j to the end
    #vertical array
    for i=j:200
        C[i] = nothing
    end
    zSet = C
    # println("zSet",zSet)
    # assigns C matrix we operated on to variabkle zSet
    T = zeros(Complex{Float64},order,2*order)
    # println("T",size(T))
    # create a new zeroes function called T
    if form =="CRFB"
        # #Find g
        # #Assume the roots are ordered, real first, then cx conj. pairs
        for i=1:order2
            # calculate g value
            g[i] = 2*(1-real(ntf_z[2*i-1+odd])) 
        end
        # # establish another matrix of 0s
        # # L1 = zeros(1,order) 
        L1 = zeros(Complex{Float64}, 1, order*2)
        #     #Form the linear matrix equation a*T*=L1
        #     #iterate
        # println("L1",size(L1))
        for i=1:order*2
            # at the i index of zSet matrix assign the value to variable z
            z = zSet[i] 
                #L1(z) = 1-1/H(z)

            y1 = ones(Complex{Float64},length(z),1)
            for i=1:length(ntf_p)
                y1 = y1.*(z-ntf_p[i])
            end
            # println("y1: ",y1)
            y2 = ones(Complex{Float64},length(z),1)
            for i=1:length(ntf_z)
                y2 = y2.*(z-ntf_z[i])
            end
            # println("y2: ",y2)
            v1=y1[1]/y2[1]
            # println("v1: ",v1)
            L1[i] = 1-v1
            println("z: ",z)
            Dfactor = (z-1)/z
            println("Dfactor: ",Dfactor)
            println("")
            product=1
            for j=order:-2:(1+odd)
                ntf_sub=ntf_z[j-1:j]
                # println("Cut vector: ",ntf_sub)
                y3 = ones(Complex{Float64},length(z),1)
                for i=1:length(ntf_sub)
                    # println("see ",z-ntf_sub[i])
                    y3 = y3.*(z-ntf_sub[i])
                    # println("y3: ",y3)
                end
                product = z/y3[1]*product
                # println("Product: ",product)
                # do some variable calculations
                # println("Product: ",product)
                # println("Dfactor: ",Dfactor)
                # println("Answer: ",Dfactor*product)
                T[j,i] = Dfactor*product
                T[j-1,i] = product
            end
            if odd==1
                T[1,i]=product/(z-1)
            end
            println("")
        end
        println("T: ",T)
        a = -real(L1/T) 
        # # get the real part of L1/T and give it a negative sign
        if isempty(stf)
            for i=1:order
                b[i]=a[i]
            end
            b[order+1] = 1
        end
    end
    println("a: ",a)
    println("b: ",b)
    println("g: ",g)
    println("c: ",c)
    return a,b,g,c
end
        # elseif form == 'CRFF'
        #      #Find g
        #      #Assume the roots are ordered, real first, then cx conj. pairs
        #     for i=1:order2
        #         g[i] = 2*(1-real(ntf_z[2*i-1+odd]))
        #     end
        #     # L1 = zeros(1,order)
        #     L1 = zeros(Float64, 1, order)
        #      #Form the linear matrix equation a*T*=L1
        #     for i=1:order*2
        #         z = zSet[i]
        #          #L1(z) = 1-1/H(z)
        #         # L1(i) = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z)
        #         L1[i] = 1-ntf_p_roots/ntf_z_roots
        #         if odd
        #             Dfactor = z-1
        #             product = 1/Dfactor
        #             T[1,i] = product
        #         else
        #             Dfactor = (z-1)/z
        #             product=1
        #         end
        #         for j=1+odd:2:order
        #             # product = z/evalRPoly(ntf_z(j:j+1),z)*product
        #             product = z/((roots(ntf_z[j:j+1]))*ones(Float64,length(z),1))*product
        #             T[j,i] = product*Dfactor
        #             T[j+1,i] = product
        #         end
        #     end
        #     a = -real(L1/T)
        #     if isempty(stf)
        #         b = [ 1 zeros(Float64,1,order-1) 1]
        #     end
            
        # else if form == 'CIFB'
        #      #Assume the roots are ordered, real first, then cx conj. pairs
        #      #Note ones which are moved significantly.
        #     if any(ntf_z->abs(real(ntf_z)-1) > 1e-3,ntf_z)
        #         # fprintf(stderr,'#s Warning: The ntf''s zeros have had their real parts set to one.\n', mfilename)
        #         println("NTF zeros had their real part set to one")
        #     end
        #     ntf_z = 1 + 1i*imag(ntf_z)
        #     for i=1:order2
        #         g[i] = imag(ntf_z[2*i-1+odd])^2
        #     end
        #     # L1 = zeros(1,order)
        #     L1 = zeros(Float64, 1, order)
        #      #Form the linear matrix equation a*T*=L1
        #     for i=1:order*2
        #         z = zSet[i]
        #          #L1(z) = 1-1/H(z)
        #         # L1(i) = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z)
        #         L1[i] = 1-ntf_p_roots/ntf_z_roots
        #         Dfactor = (z-1)
        #         product = 1
        #         for j=order:-2:(1+odd)
        #             # product = product/evalRPoly(ntf_z((j-1):j),z)
        #             product = product/((roots(ntf_z[j-1:j]))*ones(Float64,length(z),1))
        #             T[j,i] = product*Dfactor
        #             T[j-1,i] = product
        #         end
        #         if odd
        #             T[1,i] = product/(z-1)
        #         end
        #     end
        #     a = -real(L1/T)
        #     if isempty(stf)
        #         b = a
        #         b(order+1) = 1
        #     end
            
        # elseif form == 'CIFF'
        #      #Assume the roots are ordered, real first, then cx conj. pairs
        #      #Note ones which are moved significantly.
        #     if any(ntf_z->abs(real(ntf_z)-1) > 1e-3,ntf_z)
        #         println("NTF zeros had their real part set to one")
        #     end
        #     ntf_z = 1 + 1i*imag(ntf_z)
        #     for i=1:order2
        #         g[i] = imag(ntf_z[2*i-1+odd])^2
        #     end
        #     L1 = zeros(Float64, 1, order)
        #      #Form the linear matrix equation a*T*=L1
        #     for i=1:order*2
        #         z = zSet[i]
        #          #L1(z) = 1-1/H(z)
        #         # L1(i) = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z)
        #         L1[i] = 1-ntf_p_roots/ntf_z_roots
        #         Dfactor = (z-1)
        #         if odd
        #             product = 1/(z-1)
        #             T[1,i] = product
        #         else
        #             product = 1
        #         end
        #         for j=odd+1:2:order-1
        #             # product = product/evalRPoly(ntf_z(j:j+1),z)
        #             product = product/((roots(ntf_z[j:j+1]))*ones(Float64,length(z),1))
        #             T[j,i] = product*Dfactor
        #             T[j+1,i] = product
        #         end
        #     end
        #     a = -real(L1/T)
        #     if isempty(stf)
        #         b = [ 1 zeros(Float64,1,order-1) 1]
        #     end
            
    #     else if form == 'CRFBD'
    #          #Find g
    #          #Assume the roots are ordered, real first, then cx conj. pairs
    #         for i=1:order2
    #             g[i] = 2*(1-real(ntf_z[2*i-1+odd]))
    #         end
    #         L1 = zeros(Float64,1,order)
    #          #Form the linear matrix equation a*T*=L1
    #         for i=1:order*2
    #             z = zSet[i]
    #             #L1(z) = 1-1/H(z)
    #             L1(i) = 1-ntf_p_roots/ntf_z_roots
    #             Dfactor = (z-1)
    #             product=1/z
    #             for j=order:-2:(1+odd)
    #                 # product = z/evalRPoly(ntf_z((j-1):j),z)*product
    #                 product = z/((roots(ntf_z[j-1:j]))*ones(Float64,length(z),1))*product
    #                 T[j,i] = product*Dfactor
    #                 T[j-1,i] = product
    #             end
    #             if odd
    #                 T[1,i]=product*z/(z-1)
    #             end
    #         end
    #         a = -real(L1/T)
    #         if isempty(stf)
    #             b = a
    #             b(order+1) = 1
    #         end
            
    #     else if form == 'CRFFD'
    #          #Find g
    #          #Assume the roots are ordered, real first, then cx conj. pairs
    #         for i=1:order2
    #             g[i] = 2*(1-real(ntf_z[2*i-1+odd]))
    #         end
    #          #zL1 = z*(1-1/H(z))
    #         #  zL1 = zSet .* (1-1./evalTF(ntf,zSet))
    #         zL1 = zSet * (1-1\evalTF(ntf,zSet))
    #          #Form the linear matrix equation a*T*=zL1
    #         for i=1:order*2
    #             z = zSet[i]
    #             if odd
    #                 Dfactor = (z-1)/z
    #                 product = 1/Dfactor
    #                 T[1,i] = product
    #             else
    #                 Dfactor = z-1
    #                 product=1
    #             end
    #             for j=1+odd:2:order
    #                 # product = z/evalRPoly(ntf_z(j:j+1),z)*product
    #                 product = z/((roots(ntf_z[j:j+1]))*ones(Float64,length(z),1))*product
    #                 T[j,i] = product*Dfactor
    #                 T[j+1,i] = product
    #             end
    #         end
    #         a = -real(zL1/T)
    #         if isempty(stf)
    #             b = [ 1 zeros(Float64,1,order-1) 1]
    #         end
            
    #     else if form == 'PFF'
    #          #Find g
    #          #Assume the roots are ordered, real first, then cx conj. pairs
    #          # with the secondary zeros after the primary zeros
    #         for i=1:order2
    #             g[i] = 2*(1-real(ntf_z[2*i-1+odd]))
    #         end
    #          # Find the dividing line between the zeros
    #         theta0 = abs(angle(ntf_z[1]))
    #          # !! 0.5 radians is an arbitrary separation !!
    #         i = filter(ntf_z->abs(angle(ntf_z) - theta0) > 0.5,ntf_z)
    #         order_1 = i[1]-1
    #         order_2 = order-order_1
    #         if length(i) ≈ order_2
    #             println("For the PFF form, the NTF zeros must be sorted into primary and secondary zeros")
    #         end
    #         odd_1 = mod(order_1,2)
    #         odd_2 = mod(order_2,2)
    #         L1 = zeros(Float64,1,order)
    #          #Form the linear matrix equation a*T*=L1
    #         for i=1:order*2
    #             z = zSet[i]
    #              #L1(z) = 1-1/H(z)
    #             L1[i] = L1(i) = 1-ntf_p_roots/ntf_z_roots
    #             if odd_1
    #                 Dfactor = z-1
    #                 product = 1/Dfactor
    #                 T[1,i] = product
    #             else
    #                 Dfactor = (z-1)/z
    #                 product=1
    #             end
    #             for j=1+odd_1:2:order_1
    #                 # product = z/evalRPoly(ntf_z(j:j+1),z)*product
    #                 product = z/((roots(ntf_z[j:j+1]))*ones(Float64,length(z),1))*product
    #                 T[j,i] = product*Dfactor
    #                 T[j+1,i] = product
    #             end
    #             if odd_2
    #                 Dfactor = z-1
    #                 product = 1/Dfactor
    #                 T[order_1+1,i] = product
    #             else
    #                 Dfactor = (z-1)/z
    #                 product=1
    #             end
    #             for j=order_1+1+odd_2:2:order
    #                 # product = z/evalRPoly(ntf_z(j:j+1),z)*product
    #                 product = z/((roots(ntf_z[j:j+1]))*ones(Float64,length(z),1))*product
    #                 T(j,i) = product*Dfactor
    #                 T(j+1,i) = product
    #             end
    #         end
    #         a = -real(L1/T)
    #         if isempty(stf)
    #             b = [ 1 zeros(Float64,1,order_1-1) 1 zeros(Float64,1,order_2-1) 1]
    #         end
            
    #     else if 'Stratos'
    #          #Find g
    #          #Assume the roots are ordered, real first, then cx conj. pairs
    #         for i=1:order2
    #             g[i] = 2*(1-real(ntf_z[2*i-1+odd]))
    #         end
    #         L1 = zeros(Float64,1,order)
    #          #Form the linear matrix equation a*T*=L1
    #         for i=1:order*2
    #             z = zSet[i]
    #              #L1(z) = 1-1/H(z)
    #             L1[i] = L1[i] = L1(i) = 1-ntf_p_roots/ntf_z_roots
    #             Dfactor = (z-1)
    #             if ~odd
    #                 product = 1/(z-1)
    #                 T[1,i] = product
    #             else
    #                 product = 1
    #             end
    #             for j=odd+1:2:order-1
    #                 # product = product/evalRPoly(ntf_z(j:j+1),z)
    #                 product = product/((roots(ntf_z[j:j+1]))*ones(Float64,length(z),1))
    #                 T[j,i] = product*Dfactor
    #                 T[j+1,i] = product
    #             end
    #         end
    #         a = -real(L1/T)
    #         if isempty(stf)
    #             b = [ 1 zeros(Float64,1,order-1) 1]
    #         end
            
    #     case 'DSFB'
    #          #Find g
    #          #Assume the roots are ordered, real first, then cx conj. pairs
    #         for i=1:(order+odd)/2 -1
    #             g[i] = tan( angle( ntf_z(2*i-1+odd) ) /2 ) ^ 2
    #         end
    #         if ~odd
    #             g[order/2] = 2*(1 - cos( angle( ntf_z(order) ) ))
    #         end
    #         L1 = zeros(Float64,1,order)
    #          # Form the linear matrix equation a*T = L1
    #         for i=1:order*2
    #             z = zSet[i]
    #             #L1(z) = 1-1/H(z)
    #             L1[i] = 1-ntf_p_roots/ntf_z_roots
    #             if odd
    #                 product = 1/(z-1)
    #                 T[order,i] = product
    #             else
    #                 Dfactor = (z-1)^2 + g[order/2]*z
    #                 T[order,i] = (z-1)/Dfactor
    #                 product = (z+1)/Dfactor
    #                 T[order-1,i] = product
    #             end
    #             for j=order+odd-2:-2:1
    #                 Dfactor = (z-1)^2 + g[j/2]*(z+1)^2
    #                 product = product/Dfactor
    #                 T[j,i] = product*(z^2-1)
    #                 product = product*(z+1)^2
    #                 T[j-1,i] = product
    #             end
    #         end
    #         a = -real(L1/T)
    #         if isempty(stf)
    #             b = a
    #             b[order+1] = 1
    #         end
    #     else
    #         println("unsupported form")
    #     end      
    # end

    # if ~isempty(stf)
    #      # Compute the TF from each feed-in to the output
    #      # and solve for coefficients which yield the best match
    #      # THIS CODE IS NOT OPTIMAL, in terms of computational efficiency.
    #     # stfList = zeros(Float64,1,order+1)
    #     stfList = fill(-99, (1,order+1) )
    #     # create a M,N matrix of empty values
    #     for i = 1:order+1
    #         bi = zeros(Float64,1,order+1)
    #         bi[i]=1
    #         # assign 1 at designated parts
    #         ABCD = stuffABCD(a,g,bi,c,form) 
    #         # compute ABCD matrix of specifed structure
    #         # need to implement stuffABCD or try to find existing API
    #         # form is the input specified from the input
    #         # if strcmp(form(3:4),'FF')
    #         if SubString(form, 3:4)=="FF"
    #             # if forward feedback, 
    #             ABCD[1,order+2] = -1
    #             # see what the struction of ABCD look like and how the to operate on top of it
    #              # b1 is used to set the feedback
    #         end
    #         stfList[i] = calculateTF(ABCD) 
    #         # find caluclate TF funtion most likely have to implement
    #         # stfList is the empty matrix allocated, find out [ ] and what putting junk does
    #         # its doing some sort of assignment find out how it works
    #     end
    #      # Build the matrix equation b A = x and solve it
    #     A = zeros(Float64,order+1,length(zSet)) 
    #     # build A
    #     for i = 1:order+1
    #         A[i,:] = evalTF(stfList{i},zSet)
    #         evalTFTF=stfList(zSet,false)[:,1,1]
    #     end
    #     x = evalTF(stf,zSet)
    #     x = tr(x)
    #     # some Ax=b syntax and compute B
    #     b = real(x/A)
    #     # funciton returns respective a,g,b,c
#evalTF=H(z,false)[:,1,1]