function ds_optzeros(n , opt_odd)

    if n == 1
        optZeros = 0
    elseif n == 2
        if opt_odd
            optZeros = [sqrt(1/3)]
        else
            optZeros = [0]
        end
    elseif n == 3
        optZeros = [sqrt(3/5) 0]
    elseif n == 4
        if opt_odd
            discr = sqrt(9.0 / 49 - 3.0 / 35)
            tmp = 3.0 / 7
            optZeros = sqrt.([tmp+discr tmp-discr])
        else
            optZeros = [0 sqrt(5/7)]
        end
    elseif n == 5
        discr = sqrt(25.0 / 81 - 5.0 / 21)
        tmp = 5.0 / 9
        optZeros = sqrt.([tmp+discr tmp-discr 0])
    elseif n == 6
        if opt_odd
            optZeros = [ 0.23862059 0.66120988 0.9324696 ]
        else
            discr = sqrt(56.0) / 33
            tmp = 7.0/11;
            optZeros = sqrt.([0 tmp+discr tmp-discr])
        end
    elseif n == 7
        optZeros = [0 0.40584371 0.74153078 0.94910785 ]
    elseif n == 8
        if opt_odd
            optZeros = [ 0.18343709 0.52553345 0.79666684 0.96028993 ]
        else
            optZeros = [ 0 0.50563161 0.79017286 0.95914731 ]
        end
    elseif n == 9
        optZeros = [ 0 0.32425101 0.61337056 0.83603082 0.9681602 ]
    elseif n == 10
        if opt_odd
            optZeros = [ 0.1834370913 0.5255334458 0.7966668433 0.9602899327]
        else
            optZeros = [0 0.41572267 0.67208682 0.86238894 0.97342121 ]
        end
    elseif n == 11
        optZeros = [0 0.26953955 0.51909468 0.73015137 0.88706238 0.97822864]
    elseif n == 12
        if opt_odd
            optZeros = [0.12523875 0.36783403 0.58731921 0.7699033 0.90411753 0.9815607]
        else
            optZeros = [0 0.35222363 0.58006251 0.76647993 0.90281326 0.98132047]
        end
    elseif n == 13
        optZeros = [ 0 0.23045331 0.44849063 0.64234828 0.8015776 0.91759824 0.98418306 ]
    elseif n == 14
        if opt_odd
            optZeros = [ 0.10806212 0.31911586 0.51525046 0.68729392 0.82720185 0.92843513 0.98628389 ]
        else
            optZeros = [ 0 0.30524384 0.50836649 0.6836066 0.82537239 0.92772336 0.98615167 ]
        end
    else
        println("Optimized zeros for n>14 are not available.");
        return;
    end


    # Sort the zeros and replicate them.
    z = sort(optZeros[1,:])
    optZeros = zeros(n)
    m=1
    if isodd(n)
        optZeros[1] = z[1]
        z = z[2:length(z)]
        m = m + 1
    end

    for i = 1:length(z)
        optZeros[m]   =  z[i]
        optZeros[m+1] = -z[i]
        m = m + 2
    end

    return optZeros
end