using ControlSystems
using Statistics

"""
    sv, sx, sigma_se, max_sx, max_sy = simulateMS(v, M=16, mtf=zpk(1,0,1,1), d=0, dw=[1..], sx0=[0..])

Simulate the vector-based Mismatch-Shaping for a multi-element DAC.
Inputs:
- `v`:        A vector of the digital input values. `v` is in (-`M`:2:`M`) if `dw`=[1..]
              otherwise `v` is in [-sum(`dw`),sum(`dw`)].
- `M`:        The number of elements.
- `mtf`:	  The mismatch-shaping transfer function, given in zero-pole form.
- `d`:        Dither uniformly distributed in [-`d`,`d`] is added to sy.
- `dw`:       A vector of dac element weights
- `sx0`:      A matrix whose columns are the initial states of the ESL.
Outputs:
- `sv`:       An MxN matrix whose columns are the selection vectors.
- `sx`:       An orderxM matrix containing the final state of the ESL.
- `sigma_se`: The rms value of the selection error.
- `max_sx`:   The maximum absolute value of the state for all modulators.
- `max_sy`:   The maximum absolute value of the input to the VQ.

For a brief description of the theory of mismatch-shaping DACs, see
R. Schreier and B. Zhang "Noise-shaped multibit D/A convertor employing
unit elements," Electronics Letters, vol. 31, no. 20, pp. 1712-1713,
Sept. 28 1995.
"""
function simulateMS(v, M=16, mtf=zpk([1], [0], 1, 1), d=0, dw=[], sx0=[])
    order = length(mtf.matrix[1].p)
    if isempty(sx0)
        sx0 = zeros(order, M)
    end
    if isempty(dw)
        dw = ones(M)
    end

    # B/A = MTF-1
    mtf = tf(mtf)
    num = mtf.matrix[1].num[end:-1:0]
    den = mtf.matrix[1].den[end:-1:0]
    A = real(-den[2:order+1])
    B = real(num[2:order+1] + A)

    N = length(v)
    sv = zeros(Int, M, N)

    sx = sx0'
    max_sx = maximum(abs.(sx))
    max_sy = 0
    sum_se2 = 0

    for i = 1:N
        # Compute the sy vector.
        sy = sx * B
        # Normalize sy for a minimum value of zero.
        sy = sy .- minimum(sy)
        # Pick the elements that have the largest desired_usage (sy) values.
        sv[:,i] = selectElement(v[i], sy + d * (2*rand(M) .- 1), dw)
        # Compute the selection error.
        se = sv[:,i] - sy
        # Compute the new sx matrix
        sxn = sx * A + se
        sx = [sxn sx[:,1:order-1]]
        # Keep track of some statistics.
        sum_se2 = sum_se2 + sum((se .- mean(se)).^2)
        max_sx = max(max_sx, maximum(abs.(sxn)))
        max_sy = max(max_sy, maximum(abs.(sy)))
    end
    sigma_se = sqrt(sum_se2 / (M*N))

    return sv, sx', sigma_se, max_sx, max_sy
end
