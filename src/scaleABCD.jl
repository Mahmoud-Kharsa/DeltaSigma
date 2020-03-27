using LinearAlgebra
using Printf
using Random

"""
    ABCDs, umax, S = scaleABCD(ABCD, nlev=2, f=0, xlim=1, ymax=nlev+5, umax, N=100000, N0=10)

Scale the loop filter of a general delta-sigma modulator for dynamic range.

Inputs:
- `ABCD`:   The state-space description of the loop filter.
- `nlev`:   The number of levels in the quantizer.
- `xlim`:   A vector or scalar specifying the limit for each state variable.
- `ymax`:   The stability threshold. Inputs that yield quantizer inputs above `ymax`
            are considered to be beyond the stable range of the modulator.
- `umax`:   The maximum allowable input amplitude. `umax` is calculated if it
            is not supplied.
Outputs:
- `ABCDs`:  The state-space description of the scaled loop filter.
- `S`:      The diagonal scaling matrix, `ABCDs` = [`S`A`S`⁻¹ `S`B; C`S`⁻¹ D]; xs = `S`x;
"""
function scaleABCD(ABCD, nlev=2, f=0, xlim=1, ymax=NaN, umax=NaN, N_sim=100000, N0=10)
    if isnan(ymax)
        ymax = nlev + 5
    end
    order = size(ABCD, 1) - 1
    if size(xlim) == ()
        xlim = xlim * ones(order)   # Convert scalar xlim to a vector
    end

    # Make this function repeatable
    rng = MersenneTwister(0)

    # Envelope for smooth start-up
    raised_cosine = 0.5 * (1 .- cos.(pi/N0 * (0:N0-1)))
    if isnan(umax)
        # Simulate the modulator with DC or sine wave inputs to detect its
        # stable input range.
        # First get a rough estimate of umax.
        ulist = (0.1:0.1:1.0) * (nlev-1)
        umax = nlev - 1
        N = 1000
        u0 = [exp.(2im*pi*f*(-N0:-1)) .* raised_cosine; exp.(2im*pi*f*(0:N-1))]
        u0 = u0 .+ 0.01 * randn(rng, N+N0, 2) * [1, 1im]
        u0 = real(u0)

        for u = ulist
            _, _, _, y = simulateDSM(u*u0, ABCD, nlev)
            if maximum(abs.(y)) > ymax
                umax = u    # umax is the smallest input found which causes 'instability'
                break
            end
        end

        if umax == ulist[1]
            error(@sprintf("Modulator is unstable even with an input amplitude of %.1f.", umax))
        end
    end

    # More detailed simulation
    N = N_sim
    u0 = [exp.(2im*pi*f*(-N0:-1)) .* raised_cosine; exp.(2im*pi*f*(0:N-1))]
    u0 = u0 .+ 0.01 * randn(rng, N+N0, 2) * [1, 1im]
    u0 = real(u0)
    maxima = zeros(order) .- 1
    ulist = range(0.7*umax, umax, length=10)
    for u = ulist
        _, _, xmax, y = simulateDSM(u*u0, ABCD, nlev)
        if maximum(abs.(y)) > ymax
            break
        end
        # We need to do this at least once.
        umax = u    # umax is the largest input which keeps |y| < ymax
        maxima, = findmax([maxima xmax], dims=2)
    end

    # Scale the modulator so that all states are at most xlim.
    scale = maxima[:] ./ xlim
    S = diagm(1 ./ scale)
    Sinv = diagm(scale)
    A, B, C, D = partitionABCD(ABCD)
    ABCDs = [S*A*Sinv S*B; C*Sinv D]

    return ABCDs, umax, S
end
