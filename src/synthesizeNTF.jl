include("ds_optzeros.jl")
using ControlSystems

function synthesizeNTF(order, osr, opt, H_inf, f0)
    if f0 != 0
        #println("No band pass yet!")
        return
    end

    if opt == 3
        #println("No opt 3 yet!")
        return
    end 

    #Find Zeros
    dw = pi/osr
    if opt == 0
        z = zeros(order)
    else
        z = dw * ds_optzeros(order, isodd(opt))
        if isempty(z)
            return
        end
    end
    z = exp.(im*z)

    zp = z[angle.(z) .>  0]
    x0 = ( angle.(zp) .- 2 * pi * f0 ) * osr / pi
    p = zeros(order)

    ntf = zpk(z, p, 1, 1)

    Hinf_itn_limit = 100

    #Determin Poles
    opt_iteration = 5
    delta_x = 0
    fprev = 0
    while opt_iteration > 0
        ftol = 1e-10
        z_inf = -1
        
        HinfLimit = 2 ^ order
        if H_inf >= HinfLimit
            println("%s warning: Unable to achieve specified Hinf.")
            println("Setting all NTF poles to zero.")
            p = zeros(order)
        else
            x = 0.3 ^ (order-1)
            converged = 0
            for itn=1:Hinf_itn_limit

                #println("********************************************")
                me2 = -0.5 * (x^(2.0/order))
                w = (2 * collect(1:order) .- 1 ) * pi/order
                #println("w: ", w)

                mb2 = 1 .+ me2 * exp.(im * w)
                
                #println("mb2: ", mb2)
                
                for i=1:length(mb2)
                    temp1 = 0 + 0im
                    if abs(real(mb2[i])) >= 1.0e-16
                        temp1 += real(mb2[i])
                    end
                    if abs(imag(mb2[i])) >= 1.0e-16
                        temp1 += imag(mb2[i]) * im
                    end
                    mb2[i] = temp1;
                end

                #println("mb2_N: ", mb2)

                p = mb2 - sqrt.(mb2 .^ 2 .- 1)
                out = findall(x -> abs(x) > 1, p)
                p[out] = 1.0 ./ p[out]

                p = sort(p, by = x -> (abs(x), imag(x)), rev=true)

                #println("p: ", p)

                ntf = zpk(z, p, 1, 1)

                f = real(ntf(z_inf)[1, 1]) - H_inf

                

                if itn == 1
                    delta_x = -f/100;
                else
                    delta_x = -f*delta_x/(f-fprev);
                end
                
                
                xplus = x+delta_x;
                if xplus>0
                    x = xplus;
                else
                    x = x*0.1;
                end

                fprev = f;
                if (abs(f) < ftol) | (abs(delta_x ) < 1e-10)
                    converged = 1;
                    break;
                end

                if x>1e6
                    println(itn)
                    println("warning: Unable to achieve specified Hinf.")
                    println("Setting all NTF poles to zero")
                    P = zeros(order)
                    break
                end

                if itn == Hinf_itn_limit
                    println("warning: Danger! Iteration limit exceeded.")
                end
            end
        end
        if opt < 3   # Do not optimize the zeros
            opt_iteration = 0
        end
    end 
   
    ntf = zpk(z, p, 1, 1)
    #println("Passed1")
    #println(ntf);
    #println("Passed2")
    return ntf
end