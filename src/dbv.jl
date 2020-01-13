function dbv(x)
    if length(x) == 1
        if x == 0
            return -Inf
        else
            return 20 * log10(abs(x))
        end
    end

    y = -Inf * ones(size(x))
    if isempty(x)
        return []
    end
    nonzero = (x .!= 0)
    y[nonzero] = 20 * log10.(abs.(x[nonzero]))
    return y
end
