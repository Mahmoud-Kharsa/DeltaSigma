using Polynomials

"""
    a, g, b, c = realizeNTF(ntf, form="CRFB", stf=nothing)

Convert a noise transfer function into coefficients for the desired structure.

Supported structures are:
- `CRFB`:    Cascade of resonators, feedback form.
- `CRFF`:    Cascade of resonators, feedforward form.
- `CIFB`:    Cascade of integrators, feedback form.
- `CIFF`:    Cascade of integrators, feedforward form.
- `CRFBD`:   CRFB with delaying quantizer.
- `CRFFD`:   CRFF with delaying quantizer.
- `PFF`:     Parallel feed-forward.
- `Stratos`: A CIFF-like structure with non-delaying resonator feedbacks,
             contributed in 2007 by Jeff Gealow
- `DSFB`:    Cascade of double-sampled integrators, feedback form.

The order of the `ntf` zeros must be (real, complex conj. pairs).
The order of the zeros is used when mapping the `ntf` onto the chosen topology.

`stf` is a TransferFunction object.
"""
function realizeNTF(ntf, form="CRFB", stf=nothing)
    # The basic idea is to equate the loop filter at a set of
    # points in the z-plane to L1 = 1-1/ntf at those points.

    # Code common to all functions
    ntf_z = ntf.matrix[1].z
    ntf_p = ntf.matrix[1].p
    order = length(ntf_p)
    order2 = order ÷ 2
    odd = order - 2*order2

    a = zeros(order)
    g = zeros(order2)
    b = zeros(order+1)
    c = ones(order)

    # Choose a set of points in the z-plane at which to try to make L1 = 1-1/H
    # Author of MATLAB code this is based on did't know how to do this properly,
    # but the code below seems to work for order <= 11
    N = 200
    min_distance = 0.09
    C = zeros(ComplexF64, N)
    j = 1
    for i = 1:N
        z = 1.1 * exp(2im*pi*i/N)
        if all(x -> abs(x - z) > min_distance, ntf_z)
            C[j] = z
            j = j + 1
        end
    end
    zSet = C[1:j-1]
    T = zeros(ComplexF64, order, order*2)

    if form =="CRFB"
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        for i = 1:order2
            g[i] = 2 * (1 - real(ntf_z[2*i-1+odd]))
        end

        # Form the linear matrix equation a*T = L1
        L1 = zeros(ComplexF64, 1, order*2)
        for i = 1:order*2
            z = zSet[i]

            # L1(z) = 1-1/H(z)
            L1[i] = 1 - poly(ntf_p)(z) / poly(ntf_z)(z)

            Dfactor = (z - 1) / z
            product = 1
            for j = order:-2:(odd+1)
                product = z / poly(ntf_z[j-1:j])(z) * product
                T[j,i] = Dfactor * product
                T[j-1,i] = product
            end
            if odd == 1
                T[1,i] = product / (z - 1)
            end
        end

        a = -real(L1 / T)[:]
        if isnothing(stf)
            b[1:order] = a
            b[order+1] = 1
        end

    elseif form == "CRFF"
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        for i = 1:order2
            g[i] = 2 * (1 - real(ntf_z[2*i-1+odd]))
        end

        # Form the linear matrix equation a*T = L1
        L1 = zeros(ComplexF64, 1, order*2)
        for i = 1:order*2
            z = zSet[i]

            # L1(z) = 1-1/H(z)
            L1[i] = 1 - poly(ntf_p)(z) / poly(ntf_z)(z)

            if odd == 1
                Dfactor = z - 1
                product = 1 / Dfactor
                T[1,i] = product
            else
                Dfactor = (z - 1) / z
                product = 1
            end
            for j = (odd+1):2:order
                product = z / poly(ntf_z[j:j+1])(z) * product
                T[j,i] = product * Dfactor
                T[j+1,i] = product
            end
        end

        a = -real(L1 / T)[:]
        if isnothing(stf)
            b = [1; zeros(order-1); 1]
        end

    elseif form == "CIFB"
        # Assume the roots are ordered, real first, then cx conj. pairs
        # Note ones which are moved significantly.
        if any(ntf_z -> abs(real(ntf_z) - 1) > 1e-3, ntf_z)
            @warn "The ntf's zeros have had their real parts set to one."
        end
        ntf_z = 1 .+ 1im*imag(ntf_z)

        for i = 1:order2
            g[i] = imag(ntf_z[2*i-1+odd])^2
        end

        # Form the linear matrix equation a*T = L1
        L1 = zeros(ComplexF64, 1, order*2)
        for i = 1:order*2
            z = zSet[i]

            # L1(z) = 1-1/H(z)
            L1[i] = 1 - poly(ntf_p)(z) / poly(ntf_z)(z)

            Dfactor = z - 1
            product = 1
            for j = order:-2:(odd+1)
                product = product / poly(ntf_z[j-1:j])(z)
                T[j,i] = product * Dfactor
                T[j-1,i] = product
            end
            if odd == 1
                T[1,i] = product / (z - 1)
            end
        end

        a = -real(L1 / T)[:]
        if isnothing(stf)
            b[1:order] = a
            b[order+1] = 1
        end

    elseif form == "CIFF"
        # Assume the roots are ordered, real first, then cx conj. pairs
        # Note ones which are moved significantly.
        if any(ntf_z -> abs(real(ntf_z) - 1) > 1e-3, ntf_z)
            @warn "The ntf''s zeros have had their real parts set to one."
        end
        ntf_z = 1 .+ 1im*imag(ntf_z)

        for i = 1:order2
            g[i] = imag(ntf_z[2*i-1+odd])^2
        end

        # Form the linear matrix equation a*T = L1
        L1 = zeros(ComplexF64, 1, order*2)
        for i = 1:order*2
            z = zSet[i]

            # L1(z) = 1-1/H(z)
            L1[i] = 1 - poly(ntf_p)(z) / poly(ntf_z)(z)

            Dfactor = z - 1
            if odd == 1
                product = 1 / (z - 1)
                T[1,i] = product
            else
                product = 1
            end
            for  j= (odd+1):2:(order-1)
                product = product / poly(ntf_z[j:j+1])(z)
                T[j,i] = product * Dfactor
                T[j+1,i] = product
            end
        end

        a = -real(L1 / T)[:]
        if isnothing(stf)
            b = [1; zeros(order-1); 1]
        end

    elseif form == "CRFBD"
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        for i = 1:order2
            g[i] = 2 * (1 - real(ntf_z[2*i-1+odd]))
        end

        # Form the linear matrix equation a*T = L1
        L1 = zeros(ComplexF64, 1, order*2)
        for i = 1:order*2
            z = zSet[i]

            # L1(z) = 1-1/H(z)
            L1[i] = 1 - poly(ntf_p)(z) / poly(ntf_z)(z)

            Dfactor = z - 1
            product = 1 / z
            for j = order:-2:(odd+1)
                product = z / poly(ntf_z[(j-1):j])(z) * product
                T[j,i] = product * Dfactor
                T[j-1,i] = product
            end
            if odd == 1
                T[1,i] = product * z / (z - 1)
            end
        end

        a = -real(L1 / T)[:]
        if isnothing(stf)
            b[1:order] = a
            b[order+1] = 1
        end

    elseif form == "CRFFD"
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        for i = 1:order2
            g[i] = 2 * (1 - real(ntf_z[2*i-1+odd]))
        end

        # Form the linear matrix equation a*T = zL1
        zL1 = zeros(ComplexF64, 1, order*2)
        for i = 1:order*2
            z = zSet[i]

            # zL1 = z*(1-1/H(z))
            zL1[i] = z * (1 - poly(ntf_p)(z) / poly(ntf_z)(z))

            if odd == 1
                Dfactor = (z - 1) / z
                product = 1 / Dfactor
                T[1,i] = product
            else
                Dfactor = z - 1
                product = 1
            end
            for j = (odd+1):2:order
                product = z / poly(ntf_z[j:j+1])(z) * product
                T[j,i] = product * Dfactor
                T[j+1,i] = product
            end
        end

        a = -real(zL1 / T)[:]
        if isnothing(stf)
            b = [1; zeros(order-1); 1]
        end

    elseif form == "PFF"
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        # with the secondary zeros after the primary zeros
        for i = 1:order2
            g[i] = 2 * (1 - real(ntf_z[2*i-1+odd]))
        end

        # Find the dividing line between the zeros
        theta0 = abs(angle(ntf_z[1]))
        # 0.5 radians is an arbitrary separation
        i = findall(x -> abs(abs(angle(x))-theta0) > 0.5, ntf_z)
        order_1 = i[1] - 1
        order_2 = order - order_1
        if length(i) != order_2
            error("For the PFF form, the NTF zeros must be sorted into primary and secondary zeros")
        end
        odd_1 = mod(order_1, 2)
        odd_2 = mod(order_2, 2)

        # Form the linear matrix equation a*T = L1
        L1 = zeros(ComplexF64, 1, order*2)
        for i = 1:order*2
            z = zSet[i]

            # L(z) = 1-1/H(z)
            L1[i] = 1 - poly(ntf_p)(z) / poly(ntf_z)(z)

            if odd_1 == 1
                Dfactor = z - 1
                product = 1 / Dfactor
                T[1,i] = product
            else
                Dfactor = (z - 1) / z
                product = 1
            end
            for j = (odd_1+1):2:order_1
                product = z / poly(ntf_z[j:j+1])(z) * product
                T[j,i] = product * Dfactor
                T[j+1,i] = product
            end

            if odd_2 == 1
                Dfactor = z - 1
                product = 1 / Dfactor
                T[order_1+1,i] = product
            else
                Dfactor = (z - 1) / z
                product = 1
            end
            for j = (order_1+1+odd_2):2:order
                product = z / poly(ntf_z[j:j+1])(z) * product
                T[j,i] = product * Dfactor
                T[j+1,i] = product
            end
        end

        a = -real(L1 / T)[:]
        if isnothing(stf)
            b = [1; zeros(order_1-1); 1; zeros(order_2-1); 1]
        end

    elseif form == "Stratos"
        # code copied from case 'CRFF':
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        for i = 1:order2
            g[i] = 2 * (1 - real(ntf_z[2*i-1+odd]))
        end

        # code copied from case 'CIFF':
        # Form the linear matrix equation a*T = L1
        L1 = zeros(ComplexF64, 1, order*2)
        for i = 1:order*2
            z = zSet[i]

            # L1(z) = 1-1/H(z)
            L1[i] = 1 - poly(ntf_p)(z) / poly(ntf_z)(z)

            Dfactor = z - 1
            if odd == 1
                product = 1 / (z - 1)
                T[1,i] = product
            else
                product = 1
            end
            for  j= (odd+1):2:(order-1)
                product = product / poly(ntf_z[j:j+1])(z)
                T[j,i] = product * Dfactor
                T[j+1,i] = product
            end
        end

        a = -real(L1 / T)[:]
        if isnothing(stf)
            b = [1; zeros(order-1); 1]
        end

    elseif form == "DSFB"
        # Find g
        # Assume the roots are ordered, real first, then cx conj. pairs
        for i = 1:((order+odd)÷2-1)
            g[i] = tan(angle(ntf_z[2*i-1+odd]) / 2) ^ 2
        end
        if odd == 0
            g[order÷2] = 2 * (1 - cos(angle(ntf_z[order])))
        end

        # Form the linear matrix equation a*T = L1
        L1 = zeros(ComplexF64, 1, order*2)
        for i = 1:order*2
            z = zSet[i]

            # L(z) = 1-1/H(z)
            L1[i] = 1 - poly(ntf_p)(z) / poly(ntf_z)(z)

            if odd == 1
                product = 1 / (z - 1)
                T[order,i] = product
            else
                Dfactor = (z - 1)^2 + g[order÷2]*z
                T[order,i] = (z - 1) / Dfactor
                product = (z + 1) / Dfactor
                T[order-1,i] = product
            end
            for j = (order+odd-2):-2:1
                Dfactor = (z - 1)^2 + g[j÷2]*(z + 1)^2
                product = product / Dfactor
                T[j,i] = product * (z^2 - 1)
                product = product * (z + 1)^2
                T[j-1,i] = product
            end
        end

        a = -real(L1 / T)[:]
        if isnothing(stf)
            b[1:order] = a
            b[order+1] = 1
        end

    else
        error("Unsupported form $form")
    end

    if !isnothing(stf)
        # Compute the TF from each feed-in to the output
        # and solve for coefficients which yield the best match
        # THIS CODE IS NOT OPTIMAL, in terms of computational efficiency.
        stfList = Array{TransferFunction}(undef, order+1)
        for i = 1:order+1
            bi = zeros(order+1)
            bi[i] = 1
            ABCD = stuffABCD(a, g, bi, c, form)
            if occursin("FF", form)
                ABCD[1,order+2] = -1    # b1 is used to set the feedback
            end
            _, stfList[i] = calculateTF(ABCD)
        end
        # Build the matrix equation b*A = x and solve it.
        A = zeros(ComplexF64, order+1, length(zSet))
        for i = 1:order+1
            A[i,:] = stfList[i](zSet, false)[:,1,1]
        end
        x = stf(zSet, false)[:,1,1]
        x = transpose(x)
        b = real(x / A)[:]
    end

    return a, g, b, c
end
