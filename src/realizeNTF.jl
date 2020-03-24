using Pkg
Pkg.add("Roots")
using Roots
Pkg.add("Polynomials")
using Polynomials
include("evalRPoly.jl")

function realizeNTF(ntf,form,stf)
    parameters = ["ntf","form","stf"]
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