function undbv(x)
    if length(x) == 1
        return 10 ^ (x/20)
    else
        return 10 .^ (x/20)
    end
end
