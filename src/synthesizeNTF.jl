function synthesizeNTF(order = 3, osr = 64.0, opt = 0, H_inf = 1.5, f0 = 0.0)
    if opt == 3
        println("No opt 3 yet!")
        return
    end 

    # Find Zeros
    if f0 != 0 #Bandpass	
        order = order ÷ 2
        dw = pi / (2 * osr)
    else
        dw = pi / osr
    end
    
    if length(opt) == 1
        if opt == 0
            z = zeros(order)
        else
            z = dw * ds_optzeros(order, isodd(opt))
            if isempty(z)
                return
            end
            if f0 != 0
                order = order * 2
                z = z .+ 2 * pi * f0
                ztmp = [ z' ; -z' ]
                z = ztmp[:]
            end
        end
        z = exp.(im * z)
    else
        z = opt[:]
    end

    zp = z[angle.(z) .>  0]
    x0 = ( angle.(zp) .- 2 * pi * f0 ) * osr / pi
    p = zeros(order)

    ntf = zpk(z, p, 1, 1)

    Hinf_itn_limit = 100

    # Determin Poles
    opt_iteration = 5
    delta_x = 0
    fprev = 0
    while opt_iteration > 0
        ftol = 1e-10
        if f0 > 0.25
            z_inf = 1
        else
            z_inf = -1
        end
        
        if f0 == 0 # Lowpass
            HinfLimit = 2^order
            if H_inf >= HinfLimit
                println("%s warning: Unable to achieve specified Hinf.")
                println("Setting all NTF poles to zero.")
                p = zeros(order)
            else
                x = 0.3^(order - 1)
                converged = 0
                for itn = 1:Hinf_itn_limit

                    # println("********************************************")
                    me2 = -0.5 * (x^(2.0 / order))
                    w = (2 * collect(1:order) .- 1 ) * pi / order
                    # println("w: ", w)

                    mb2 = 1 .+ me2 * exp.(im * w)
                    
                    # println("mb2: ", mb2)

                    # println("mb2_N: ", mb2)

                    p = mb2 - sqrt.(mb2.^2 .- 1)
                    out = findall(x->abs(x) > 1, p)
                    p[out] = 1.0 ./ p[out]

                    #println("u: ", p)

                    p = sort(p, by = x->(abs(x), imag(x)), rev = true)

                    ntf = createZPK(z, p, 1, 1)

                    f = real(ntf(z_inf)[1, 1]) - H_inf
                    #println("f: ", f)
                    

                    if itn == 1
                        delta_x = -f / 100
                    else
                        delta_x = -f * delta_x / (f - fprev)
                    end
                    
                    
                    xplus = x + delta_x
                    if xplus > 0
                        x = xplus
                    else
                        x = x * 0.1
                    end

                    fprev = f
                    if (abs(f) < ftol) | (abs(delta_x) < 1e-10)
                        converged = 1
                        break
                    end

                    if x > 1e6
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
        else #bandpass
            x = 0.3 ^ (order ÷ 2 - 1)   #starting guess (not very good for f0~0)
            c2pif0 = cos(2 * pi * f0)

            for itn=1:Hinf_itn_limit
                e2 = 0.5 * x ^ (2/order)
                w = (2 * collect(1:order) .- 1) * pi / order
                mb2 = c2pif0 .+ e2 * exp.(im * w)
                p = mb2 - sqrt.(mb2 .^ 2 .- 1)

                # reflect poles to be inside the unit circle.
                out = findall(x->abs(x) > 1, p)
                p[out] = 1.0 ./ p[out]

                p = sort(p, by = x->(abs(x), imag(x)), rev = true)

                ntf = createZPK(z, p, 1, 1)

                f = real(ntf(z_inf)[1, 1]) - H_inf

                println("f: ", f)

                # 	[x f]
                if itn==1
                    delta_x = -f / 100
                else
                    delta_x = -f * delta_x / (f - fprev)
                end
                
                xplus = x + delta_x
                if xplus > 0
                    x = xplus
                else
                    x = x * 0.1
                end
                fprev = f
                if (abs(f) < ftol) | (abs(delta_x) < 1e-10)
                    break
                end

                if x>1e6
                    println("warning: Unable to achieve specified Hinf.\n")
                    println("Setting all NTF poles to zero.\n")
                    P = zeros(order)
                    break
                end

                if itn == Hinf_itn_limit
                    println("warning: Danger! Hinf iteration limit exceeded.\n")
                end
            end
        end



        #optimization
        if opt < 3   # Do not optimize the zeros
            opt_iteration = 0
        else
            return
        end
    end 
   
    ntf = createZPK(z, p, 1, 1)
    return ntf
end