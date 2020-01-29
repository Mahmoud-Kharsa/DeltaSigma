function createZPK(num, den, k, Ts=0)

    #println("p: ", den)
    p_temp = copy(den)
    for i = 1:length(den)
        p_temp[i] = 0 + 0im
        p_temp[i] += round(real(den[i]); digits = 16)
        p_temp[i] += round(imag(den[i]); digits = 16) * im
    end
    #println("p_new: ", p_temp)

    i = 1
    while i <= length(den)
        #println("i: ", i)
        if imag(p_temp[i]) != 0
            if imag(p_temp[i]) != imag(p_temp[i+1]) || real(p_temp[i]) != real(p_temp[i+1])
                p_temp[i] = p_temp[i+1] - 2 * imag(p_temp[i+1]) * im
                if imag(p_temp[i]) < imag(p_temp[i+1])
                    p_temp[i] -= 2 * imag(p_temp[i]) * im
                    p_temp[i+1] -= 2 * imag(p_temp[i+1]) * im
                end
                i += 1
            end
        end
        i += 1
    end
    #println("p_new2: ", p_temp)

    ntf = zpk(num, p_temp, k, Ts)

    return ntf
    
end