"""
    nx = dsmap(u, ABCD, nlev, x, e, v)

For a DSM with input `u`, a structure `ABCD` and an `nlev`-level quantizer,
compute the (potential) vertices of the image of a convex object described
in terms of its vertices `x` and edges `e`. `v` is the assumed quantizer
output; it is computed if it is not supplied.

Basic Algorithm:
1) Compute the images of the vertices.
2) For those edges which cross a splitting hyperplane,
   compute the images of each split point.
"""
function dsmap(u, ABCD, nlev, x, e, v=nothing)
    n = size(ABCD, 1) - 1
    A = ABCD[1:n, 1:n]
    B = ABCD[1:n, n+1:n+2]
    C = ABCD[n+1, 1:n]
    D1 = ABCD[n+1, n+1] # For any DSM, D2=ABCD(n+1,n+2) must be zero.

    N = size(x, 2)

    # Compute v. The assumption that D1=0 for u-ranges is implicit in this step.
    if isnothing(v)
        y = x'*C .+ D1*u
        v = ds_quantize.(y, nlev)
    elseif size(v) == ()
        v = v * ones(Int, N)
    elseif size(v, 2) != N
        error("the supplied v argument is the wrong size.")
    end

    # 1) Compute the images of the vertices.
    B1u = B[:,1] * u
    nx = A*x + B1u*ones(1, N) + B[:,2]*v'

    # 2) For those edges which cross a (or several) splitting hyperplanes,
    #    compute the two images of each split point.
    diff = abs.(v[e[1,:]] - v[e[2,:]])
    split1 = diff .== 2	# edges split in one place only

    # Handle the split1 edges en masse.
    if any(split1)
        i1 = e[1,split1]
        i2 = e[2,split1]
        y0 = 0.5 * (v[i1] + v[i2]) # the appropriate quantizer thresholds
        k1 = (y[i2] - y0) ./ (y[i2] - y[i1])
        k2 = 1 .- k1
        psplit = ones(n)*k1' .* x[:,i1] + ones(n)*k2' .* x[:,i2]
        N = length(k1)
        images1 = A*psplit + B1u*ones(1,N) + B[:,2]*v[i1]'
        images2 = images1 + B[:,2]*(v[i2] - v[i1])'
        nx = [nx images1 images2]
    end

    return nx
end
