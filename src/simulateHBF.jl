using ControlSystems

"""
    y = simulateHBF(x, f1, f2, mode=0)

Simulate a Saramaki half-band filter.

The `f1` and `f2` vectors contain coefficients for the structure. (`f1` and `f2`
can each also be a tuple of arrays (val, csd) where the first array contains the
coeffiecients).

The `mode` flag determineswhether the input is filtered, interpolated, or
decimated according to the following table:
- `mode` = 0:   Plain filtering, no interpolation or decimation.
- `mode` = 1:   The input is interpolated.
- `mode` = 2:   The output is decimated, even samples are taken.
- `mode` = 3:   The output is decimated, odd samples are taken.
"""
function simulateHBF(x, f1=NaN, f2=NaN, mode=0)
    if f1 isa Number && isnan(f1)
        f1, f2 = exampleHBF(4)
    end
    if f1 isa Tuple # Presume that f1 is a (val, csd) tuple
        f1 = f1[1]
    end
    if f2 isa Tuple # Presume that f2 is a (val, csd) tuple
        f2 = f2[1]
    end

    n1 = length(f1)
    n2 = length(f2)
    f2imp = [f2[end:-1:1]; f2]
    if mode == 0    # Insert zeros
        f2imp = [f2imp'; zeros(1, 2*n2)]
        f2imp = f2imp[1:4*n2-1]
    end
    F2 = tf(f2imp, [1; zeros(length(f2imp)-1)], 1)

    if mode == 0    # Plain
        up, = lsim(F2, x, 1:length(x))[:]
        y = 0.5 * delay(x, 2*n2-1) + f1[1] * up
        for i = 2:n1
            up, = lsim(F2, up, 1:length(up))
            up, = lsim(F2, up, 1:length(up))
            y = f1[i] * up + delay(y, 4*n2-2)
        end

    elseif mode == 1    # Interpolating
        up = zeros(size(x))
        nz = 2*n2 - 1
        for i = n1:-1:1
            if i == 1
                up, = lsim(F2, up + f1[i] * x, 1:length(x))
                nz = n2 - 1
            else
                up, = lsim(F2, up + f1[i] * x, 1:length(x))
                up, = lsim(F2, up, 1:length(up))
            end
            x = delay(x, nz)
        end
        y = [2*up'; x']
        y = y[:]    # Interleave the upper and lower streams

    elseif mode == 2 || mode == 3   # Decimating
        x = x[1:2*(end÷2)]
        if mode == 3
            y = 0.5 * x[1:2:end]
            up = x[2:2:end]
            nz = n2 - 1
        else
            y = 0.5 * x[2:2:end]
            up = x[1:2:end]
            nz = n2
        end
        for i = 1:n1
            if i == 1
                up, = lsim(F2, up, 1:length(up))
            else
                up, = lsim(F2, up, 1:length(up))
                up, = lsim(F2, up, 1:length(up))
            end
            y = f1[i] * up + delay(y, nz)
            nz = 2 * n2 - 1
        end

    else
        throw(ArgumentError("$mode is not a valid value for the mode variable."))
    end

    return y[:]
end
