using ControlSystems

"""
    ntf = synthesizeNTF(order=3, osr=64, opt=0, H_inf=1.5, f0=0)

Synthesize a noise transfer function for a delta-sigma modulator.
- `order`: order of the modulator
- `osr`: oversampling ratio
- `opt`: flag for optimized zeros
    - 0 -> not optimized
    - 1 -> optimized
    - 2 -> optimized with at least one zero at band-center
    - [z] -> zero locations in complex form
- `H_inf`: maximum NTF gain
- `f0`: center frequency (1->fs)

`ntf` is a TransferFunction object containing the zeros and poles of the NTF.
"""
function synthesizeNTF(order=3, osr=64.0, opt=0, H_inf=1.5, f0=0.0)
    if f0 > 0.5
        error("f0 must be less than 0.5")
    end
    if f0 != 0 && f0 < (0.25/osr)
        @warn "Creating a lowpass ntf"
        f0 = 0
    end
    if f0 != 0 && isodd(order)
        error("order must be even for a bandpass modulator")
    end
    if length(opt) > 1 && length(opt) != order
        error("The opt vector must be of length $order(=order)")
    end

    # Determine the zeros
    if f0 != 0  # Bandpass design-- halve the order temporarily
        order = order ÷ 2
        dw = pi / (2*osr)
    else
        dw = pi / osr
    end

    if length(opt) == 1
        if opt == 0
            z = zeros(order)
        else
            z = dw * ds_optzeros(order, opt)
            if isempty(z)
                return
            end
            if f0 != 0  # Bandpass design -- shift and replicate the zeros.
                order = order * 2
                z = z .+ 2*pi*f0
                z = [z'; -z'][:]
            end
        end
        z = exp.(im * z)
    else
        z = opt[:]
    end

    p = zeros(order)
    ntf = zpk(z, p, 1, 1)
    itn_limit = 100

    delta_x = 0
    fprev = 0

    # Iteratively determine the poles by finding the value of the x-parameter
    # which results in the desired H_inf.
    if f0 == 0  # Lowpass design
        HinfLimit = 2^order
        if H_inf >= HinfLimit
            @warn "Unable to achieve specified Hinf. Setting all NTF poles to zero."
            p = zeros(order)
        else
            x = 0.3 ^ (order-1)   # starting guess
            for itn = 1:itn_limit
                me2 = -0.5 * x^(2/order)
                w = (2*(1:order) .- 1) * pi / order
                mb2 = 1 .+ me2 * exp.(im * w)
                p = mb2 - sqrt.(mb2.^2 .- 1)
                out = findall(x -> abs(x) > 1, p)
                p[out] = 1 ./ p[out] # reflect poles to be inside the unit circle
                p = pairComplex(p)

                ntf = zpk(z, p, 1, 1)
                f = real(ntf(-1, false)[1, 1]) - H_inf
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
                if abs(f) < 1e-10 || abs(delta_x) < 1e-10
                    break
                end

                if x > 1e6
                    @warn "Unable to achieve specified Hinf. Setting all NTF poles to zero."
                    p = zeros(order)
                    break
                end

                if itn == itn_limit
                    @warn "Danger! Iteration limit exceeded."
                end
            end
        end
    else    # Bandpass design
        x = 0.3 ^ (order/2 - 1)   # starting guess (not very good for f0~0)
        if f0 > 0.25
            z_inf = 1
        else
            z_inf = -1
        end
        c2pif0 = cos(2 * pi * f0)

        for itn = 1:itn_limit
            e2 = 0.5 * x^(2/order)
            w = (2*(1:order) .- 1) * pi / order
            mb2 = c2pif0 .+ e2 * exp.(im * w)
            p = mb2 - sqrt.(mb2.^2 .- 1)
            out = findall(x -> abs(x) > 1, p)
            p[out] = 1 ./ p[out]    # reflect poles to be inside the unit circle.
            p = pairComplex(p)

            ntf = zpk(z, p, 1, 1)
            f = real(ntf(z_inf, false)[1, 1]) - H_inf
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
            if abs(f) < 1e-10 || abs(delta_x) < 1e-10
                break
            end

            if x > 1e6
                @warn "Unable to achieve specified Hinf. Setting all NTF poles to zero."
                p = zeros(order)
                break
            end

            if itn == itn_limit
                @warn "Danger! Iteration limit exceeded."
            end
        end
    end

    ntf = zpk(z, p, 1, 1)
    return ntf
end

"""
    pairComplex(a)

A helper function to pair up and order complex conjugate pairs in `a`
so it can be passed to `zpk`(). Based off MATLAB cplxpair().
"""
function pairComplex(a)
    a = sort(a, by=x->(real(x), imag(x)), rev=true)
    rtol = 100 * eps(Float64) # tolerance
    i = 1
    while i <= length(a)
        if isapprox(a[i], real(a[i]), rtol=rtol)
            # imaginary component is zero or very small, treat as purely real
            a[i] = real(a[i])
        else
            if isapprox(a[i], a[i+1]', rtol=rtol)
                # complex conjugate pairs with some error, turn to exact pairs
                # by setting real and imaginary to midpoint of the two
                a_real = (real(a[i]) + real(a[i+1])) / 2
                a_imag = (abs(imag(a[i])) + abs(imag(a[i+1]))) / 2
                a[i] = a_real + a_imag*im   # put positive pair first
                a[i+1] = a_real - a_imag*im
            else
                error("Complex numbers can't be paired")
            end
            i += 1
        end
        i += 1
    end
    return a
end
