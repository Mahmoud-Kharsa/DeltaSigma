"""
    sv = selectElement(v, sy, dw=[1,1,...,1], tri=0)

Select elements of a multi-element DAC to minimize the selection error,
subject to the constraint that the nominal DAC output is `v`, i.e. `v` = (`sv`)ᵀ(`dw`).
Assume that the preferred usage order is that given by `dw`.

Input
- `v`:    DSM output in -M:2:M
- `sy`:   M×1 desired usage vector
- `dw`:   M×1 vector of element weights
- `tri`:  Flag indicating the elements are tri-level, i.e. `sv` = 0 is allowed. NOT SUPPORTED YET
Output
- `sv`:   M×1 selection vector. ±1 if `tri` = 0, {0,±1} if `tri` = 1.
"""
function selectElement(v, sy, dw=[], tri=0)
    M = length(sy)
    if isempty(dw) || all(dw .== 1)
        dw_is_one = true
        sdw = M
    else
        dw_is_one = false
        sdw = sum(dw)
    end

    if abs(v) > sdw
        error("|v| is too large.");
    elseif v == sdw
        return ones(M)
    elseif v == -sdw
        return -ones(M)
    end

    sv = -ones(Int, M)
    possibilities = sortperm(-sy)  # Determine the element priority

    # Speed up selection in the usual case where all element weights are one.
    # Suggested 1997-03-05 by J. A. Cherry.
    if dw_is_one
        sv[possibilities[1:(M+v)÷2]] .= 1
        return sv
    end

    # Go through sv possibilities one by one, until one which meets the
    # v = sv'*dw constraint is found.
    i = 1                    # Selection level
    pointer = ones(Int, M)   # Array of pointers to selected elements
    selected = zeros(Int, M) # Selected elements
    dw0 = [0; dw]
    while true
        while pointer[i] > M
            # backtrack
            i = i - 1
            if i == 0
                break   # failure!
            end
            pointer[i] = pointer[i] + 1;
            selected[i+1:end] .= 0
        end
        selected[i] = possibilities[pointer[i]]
        dv = 2 * sum(dw0[selected.+1]) - sdw
        if dv == v
            break       #success!
        elseif dv < v
            # Proceed to the next level of selection
            pointer[i+1] = pointer[i] + 1
            i = i + 1
        else
            # Try the next element at the current level.
            pointer[i] = pointer[i] + 1
        end
    end
    selected = selected[selected .!= 0]
    sv[selected] .= 1

    return sv
end
