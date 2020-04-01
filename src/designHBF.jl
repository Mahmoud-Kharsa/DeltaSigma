using DSP
using Polynomials
using Printf

"""
    f1, f2, info = designHBF(fp=0.2, delta=1e-5, debug=0)

Design a half-band filter which can be realized without general multipliers.
The filter is a composition of a prototype and sub-filter.

Input:
- `fp`:     The normalized cutoff frequency of the filter. Due to the symmetry
            imposed by a HBF, the stopband begins at 0.5-`fp`.
- `delta`:  The absolute value of the deviation of the frequency response from
            the ideal values of 1 in the passband and 0 in the stopband.

Output:
- `f1`,`f2`:    The coefficients of the prototype and sub-filters and their
                canonical-signed digit (csd) representation.
- `info`:       A vector containing the following data (only set when `debug`=1):
    - complexity:   The number of additions per output sample.
    - n1,n2:        The length of the `f1` and `f2` vectors.
    - sbr:          The achieved stob-band attenuation (dB).
    - phi:          The scaling factor for the F2 filter.
"""
function designHBF(fp=0.2, delta=1e-5, debug=false)
    worse = 0
    f1_saved = nothing
    f2_saved = nothing
    phi_saved = nothing

    # Try several different values for the fp1 parameter.
    # The best values are usually around .04
    # Surrender if 3 successive attempts yield progressively greater complexity.
    lowest_complexity = Inf
    prev_complexity = Inf
    for fp1 = [.03, .035, .025, .040, .020, .045, .015, .05]
        failed = false

        f1, zetap, phi = designF1(delta, fp1)
        if zetap == 1	# designF1 failed
            failed = true
            if debug
                @printf("designF1 failed at fp1=%f\n", fp1)
            end
        end

        if !failed
            f2 = designF2(fp, zetap, phi)
            n1 = length(f1[1])
            n2 = length(f2[1])
            if n2 == 0  # designF2 failed
                failed = true
                if debug
                    @printf("designF2 failed when zetap=%f, phi=%f\n", zetap, phi)
                end
            end
        end

        if !failed
            f1_csd_size = sum(map(x -> size(x, 2), f1[2]))
            f2_csd_size = sum(map(x -> size(x, 2), f2[2]))
            complexity = f1_csd_size + (2*n1 - 1)*(n2 + f2_csd_size - 1)

            if debug
                msg = @sprintf("%d adders: n1=%d, n2=%d, (fp1=%.2f, zetap=%.3f, phi=%4.2f)",
                                complexity, n1, n2, fp1, zetap, phi)
            else
                msg = ""
            end

            fresp, pbr, sbr = frespHBF([], f1, f2, phi, fp, msg)
            if pbr <= delta && sbr <= delta
                complexity = complexity + sbr
                if complexity < prev_complexity
                    worse = 0
                    if complexity < lowest_complexity
                        lowest_complexity = complexity
                        f1_saved = f1
                        f2_saved = f2
                        phi_saved = phi
                        if debug
                            @printf("%s\n", msg)
                        end
                    end
                else
                    worse = worse + 1
                    if worse > 2
                        break
                    end
                end
                prev_complexity = complexity
            end # if phr <= delta
        end
    end # for fp1

    if isinf(lowest_complexity)
        @printf("Unable to meet the design requirements.\n")
    else
        complexity = floor(lowest_complexity)
        msg = @sprintf("Final Design: %d adders", complexity)
        _, pbr, sbr = frespHBF([], f1_saved, f2_saved, phi_saved, fp, msg)
        n1 = length(f1_saved[1])
        n2 = length(f2_saved[1])
        if debug
            @printf("%s (%d,%d,%.0fdB)\n", msg, n1, n2, amp2db(sbr))
        end
        info = [complexity; n1; n2; amp2db(sbr); phi_saved]
    end

    return f1_saved, f2_saved, info
end

"""
    f1, zetap, phi = designF1(delta, fp1)

Design the F1 sub-filter of a Saramaki halfband filter.
This function is called by designHBF.

- `f1`:  A 2-tuple of arrays containing the F1 filter coefficents and their CSD
         representation.
- `phi`: The scaling factor for the F2 filter (imbedded in the `f1` coeffs).
"""
function designF1(delta, fp1)
    f1_saved = nothing
    zetap = 0
    phi = 0.0

    passband = exp.(4*pi*im*range(0, fp1, length=100))
    ok = false

    h = []
    n1 = 0
    for i = 1:2:7  # Odd values only
        n1 = i
        if n1 == 1
            h = [0.5, 0.5]
        else
            h = remez(2*n1, [0, 4*fp1], [1], Hz=2)
            if !(abs(sum(h) - 1) < 1e-3)   # remex bug! use least squares instead
                # firls is an open issue in DSP.jl (https://github.com/JuliaDSP/DSP.jl/issues/109)
                error("remez bug not handled yet, try again")
            end
        end
        fresp = abs.(Poly(h)(passband))
        if maximum(abs.(fresp .- 1)) <= delta
            ok = true
            break
        end
    end
    if !ok
        zetap = 1   # Use this as an indication that the function failed.
        return f1_saved, zetap, phi
    end

    # Transform h(n) to a chebyshev polynomial f1(n)
    # Sum(f1(i)*cos(w)^n)|i=1:n1 + Sum(h(n1+i))*cos(n*w))|i=1:n1, n = 2*i-1;
    w = pi * randn(n1)
    cos_w = cos.(w)
    A = zeros(n1, length(w))
    B = zeros(n1)
    for i = 1:n1
        n = 2*i - 1
        A[i,:] = cos_w .^ n
        B = B + h[n1+i] * cos.(n*w)
    end
    f1 = (transpose(B) / A)[:]

    phivecb = []

    # Optimize the quantized version of f1 to maximize the stopband width
    # ( = acos(zetap) )

    zetap = 1
    testPoints = [0; exp10.(range(-2, 0, length=128))] .- 1
    for nsd = 3:8
        # First try the unperturbed filter.
        f1a = f1
        f1b = f1
        for phia = 1 ./ [1; f1]
            phia = phia / 2.0^ceil(log2(abs(phia))) # keep phi in (0.5,1]
            # Try a bunch of coefficients in the current neighborhood,
            # shrinking the neighborhood once 10 successive trial values show no
            # improvement.  If 2 successive shrinkages do no good, try a higher nsd.
            count = 0
            nohelp = 0
            neighborhood = .05
            while neighborhood > 1e-5
                phivec = phia .^ (1:2:(2*n1) .- 1)
                if isempty(phivecb)
                    phivecb = phivec
                end
                f1q = bquantize(f1a .* phivec, nsd)
                F1 = evalF1(f1q[1], testPoints, phia)
                fi = findall(x -> abs(x) > delta, F1)
                zeta = -testPoints[max(fi[1]-1, 1)]
                if zeta < zetap
                    count = 0
                    nohelp = 0
                    zetap = zeta
                    f1b = f1q[1]
                    f1_saved = f1q
                    phi = phia
                    phivecb = phivec
                else
                    count = count + 1;
                end
                if count > 10
                    count = 0
                    neighborhood = neighborhood / 2
                    nohelp = nohelp + 1
                    if nohelp > 2
                        break
                    end
                end
                f1a = f1b./phivecb + neighborhood*(randn(size(f1b)) .- 0.5)
                phia = phia + neighborhood*(rand() - 0.5)
            end
            if zetap < 1    # Found a filter with adequate attn.
                break
            end
        end # for phia ...
        if zetap < 1    # Found a filter with adequate attn.
            break
        end
    end

    return f1_saved, zetap, phi
end

"""
    f2 = designF2(fp, zetap, phi)

Design the F2 sub-filter of a Saramaki halfband filter.
This function is called by designHBF.

subfilter design:
-    1 - delta2' < |F2/`phi`| < 1   for f in [0 `fp`];
-   -1 < |F2/`phi`| < -1 + delta2'  for f in [0.5-`fp`, 0.5];
-    1-delta2' = (1-delta2)/(1+delta2)
"""
function designF2(fp, zetap, phi)
    delta2 = (1 - zetap) / (1 + zetap)

    # determine the minimum order required by the filter
    passband = exp.(im * range(0, 4*pi*fp, length=100))
    nsub = 0
    for i = 3:2:17
        nsub = i
        h2 = remez(nsub+1, [0, 4*fp], [1], Hz=2)
        mag = abs.(Poly(h2)(passband))
        if maximum(abs.(mag .- 1)) < delta2
            break
        end
    end
    n2min = (nsub + 1) / 2

    # Search all n2,nsd pairs, in order of the product n2*(nsd+1)
    # allowing fp to be a variable?
    success = false
    nsdmin = 3
    nsdmax = 6
    for product = (nsdmin+1)*n2min:(nsdmax+1)*n2min
        for nsd = nsdmin:nsdmax
            n2 = product / (nsd + 1)
            if floor(n2) != n2  # Only take integer n2,nsd pairs
                break
            end
            n2 = floor(Int, n2)
            nsub = 2*n2 - 1

            h2 = remez(nsub+1, [0, 4*fp], [1], Hz=2)
            h2 = h2 / (phi * (1 + delta2))  # Adjust the coefficients.
            f2 = bquantize(h2[n2+1:nsub+1], nsd)
            h2 = (1 + delta2) * phi * [f2[1][n2:-1:1]; f2[1]]
            mag = abs.(Poly(h2)(passband))
            if maximum(abs.(mag .- 1)) < delta2
                success = true
                break
            end
        end
        if success
            break
        end
    end
    if !success
        f2 = ([], [])
        q2 = []
    end

    return f2
end
