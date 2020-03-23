using Pkg
Pkg.add("Roots")
using Roots
Pkg.add("Polynomials")
using Polynomials
include("evalRPoly.jl")

function realizeNTF(ntf,form,stf)
    parameters = ["ntf","form","stf"]
    defaults = [NaN, "CRFBD", []]
    ntf_z = [1.00000000000000 + 0.00000000000000im;0.999188974764100 + 0.0402665209555813im;0.999188974764100 - 0.0402665209555813im;0.997703673259689 + 0.0677302027468142im;0.997703673259689 - 0.0677302027468142im]
    ntf_p = [0.807058252221672 - 0.119565829763585im 0.807058252221672 + 0.119565829763585im 0.898412456856058 - 0.219273016428406im 0.898412456856058 + 0.219273016428406im 0.778301510991778 + 0.00000000000000im]
    
    order = length(ntf_p)
    order2 = order ÷ 2
    odd = order - 2*order2
    a = zeros(Float64, 1, order)
    g = zeros(Float64, 1, order2)
    b = zeros(Float64, 1, order+1)
    c = ones(Float64, 1, order)
    stf=[]
    N = 200
    min_distance = 0.09
    C = zeros(ComplexF64, N, 1)
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
    for i=j:200
        C[i] = nothing
    end
    zSet = C
    T = zeros(Complex{Float64},order,2*order)
    if form =="CRFB"
        for i=1:order2
            g[i] = 2*(1-real(ntf_z[2*i-1+odd])) 
        end
        L1 = zeros(Complex{Float64}, 1, order*2)
        for i=1:order*2
            z = zSet[i]
            L1[i] = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z)
            Dfactor = (z-1)/z
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
        end
        a = -real(L1/T) 
        if isempty(stf)
            for i=1:order
                b[i]=a[i]
            end
            b[order+1] = 1
        end
    elseif form == "CRFF"
        for i=1:order2
            g[i] = 2*(1-real(ntf_z[2*i-1+odd]))
        end
        L1 = zeros(Complex{Float64}, 1, order*2)
        for i=1:order*2
            z = zSet[i]
            L1[i] = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z)
            if odd==1
                Dfactor = z-1
                product = 1/Dfactor
                T[1,i] = product
            else
                Dfactor = (z-1)/z
                product=1
            end
            for j=1+odd:2:order
                product = z/evalRPoly(ntf_z[j:j+1],z)*product
                T[j,i] = product*Dfactor
                T[j+1,i] = product
            end
        end
        a = -real(L1/T);
        if isempty(stf)
            b = hcat(1,zeros(Int64,1,order-1),1)
        end
    elseif form == "CIFB"
        if any(ntf_z->abs(real(ntf_z)-1) > 1e-3,ntf_z)
            println("NTF zeros had their real part set to one")
        end
        for i=1:length(ntf_z)
            ntf_z[i]=1+1im*imag(ntf_z[i])
        end
        for i=1:order2
            g[i] = imag(ntf_z[2*i-1+odd])^2
        end
        L1 = zeros(Complex{Float64}, 1, order*2)
        for i=1:order*2
            z = zSet[i]
            L1[i] = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z)
            Dfactor = (z-1)
            product = 1
            for j=order:-2:(1+odd)
                product = product/evalRPoly(ntf_z[(j-1):j],z)
                T[j,i] = product*Dfactor
                T[j-1,i] = product
            end
            if odd==1
                T[1,i] = product/(z-1)
            end
        end
        a = -real(L1/T)
        if isempty(stf)
            for i=1:order
                b[i]=a[i]
            end
            b[order+1] = 1
        end
    elseif form == "CIFF"
        if any(ntf_z->abs(real(ntf_z)-1) > 1e-3,ntf_z)
            println("NTF zeros had their real part set to one")
        end
        for i=1:length(ntf_z)
            ntf_z[i]=1+1im*imag(ntf_z[i])
        end
        for i=1:order2
            g[i] = imag(ntf_z[2*i-1+odd])^2
        end
        L1 = zeros(Complex{Float64}, 1, order*2)
        for i=1:order*2
            z = zSet[i]
            L1[i] = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z)
            Dfactor = (z-1)
            if odd==1
                product = 1/(z-1)
                T[1,i] = product
            else
                product = 1
            end
            for j=odd+1:2:order-1
                product = product/evalRPoly(ntf_z[j:j+1],z)
                T[j,i] = product*Dfactor
                T[j+1,i] = product
            end
        end
        a = -real(L1/T)
        if isempty(stf)
            b = hcat(1,zeros(Int64,1,order-1),1)
        end
    elseif form == "CRFBD"
        for i=1:order2
            g[i] = 2*(1-real(ntf_z[2*i-1+odd]))
        end
        L1 = zeros(Complex{Float64},1,order*2)
        for i=1:order*2
            z = zSet[i]
            L1[i] = 1-evalRPoly(ntf_p,z)/evalRPoly(ntf_z,z)
            Dfactor = (z-1)
            product=1/z
            for j=order:-2:(1+odd)
                product = z/evalRPoly(ntf_z[(j-1):j],z)*product
                T[j,i] = product*Dfactor
                T[j-1,i] = product
            end
            if odd==1
                T[1,i]=product*z/(z-1)
            end
        end
        a = -real(L1/T)
        if isempty(stf)
            for i=1:order
                b[i]=a[i]
            end
            b[order+1] = 1
        end
    else
    end
    println("a: ",a)
    println("b: ",b)
    println("g: ",g)
    println("c: ",c)
    return a,b,g,c
end       
    # elseif form == "CRFFD"
    #     for i=1:order2
    #         g[i] = 2*(1-real(ntf_z[2*i-1+odd]))
    #     end
    #     zL1 = zSet * (1-1\evalTF(ntf,zSet))
    #     for i=1:order*2
    #         z = zSet[i]
    #         if odd
    #             Dfactor = (z-1)/z
    #             product = 1/Dfactor
    #             T[1,i] = product
    #         else
    #             Dfactor = z-1
    #             product=1
    #         end
    #         for j=1+odd:2:order
    #             # product = z/evalRPoly(ntf_z(j:j+1),z)*product
    #             product = z/((roots(ntf_z[j:j+1]))*ones(Float64,length(z),1))*product
    #             T[j,i] = product*Dfactor
    #             T[j+1,i] = product
    #         end
    #     end
    #     a = -real(zL1/T)
    #     if isempty(stf)
    #         b = [ 1 zeros(Float64,1,order-1) 1]
    #     end
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