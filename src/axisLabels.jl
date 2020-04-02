using Printf

function axisLabels(range, incr)
    range[abs.(range) .< 1e-6] .= 0

    s = Vector{String}(undef, length(range))
    for i = 1:length(range)
        s[i] = ""
    end

    for i = 1:incr:length(range)
        s[i] = @sprintf("%g", range[i])
    end

    return s
end
