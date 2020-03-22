using ControlSystems
using Polynomials

function simulateMS(v, M, mtf, d, dw, sx0)

    dw = []
    sx0 = []
    M = 16
    d = 0 

    order = length(mtf.matrix[1].p);
    if isempty(dw)
        dw = ones(M, 1);
    end
    if isempty(sx0)
        sx0 = zeros(order, M);
    end


    # B/A = MTF-1
    num = poly(mtf.matrix[1].z);
    den = poly(mtf.matrix[1].p);

    A = real(-den.a[end-1:-1:1]); 
    B = real(num.a[end-1:-1:1] + A);
    A = A'
    B = B'

    N = length(v);
    sv = zeros(M,N);
    
    sx = sx0;
    max_sx = maximum(abs.(sx));
    max_sy = 0;
    sum_se2 = 0;
    
    for i = 1:N
        # Compute the sy vector.
        sy = B*sx;
        # Normalize sy for a minimum value of zero.
        sy = sy .- minimum(sy);
        # Pick the elements that have the largest desired_usage (sy) values.
        sv[:,i] = selectElement(v[i], sy + d * (2 * rand(1,M) .- 1), dw);
        # Compute the selection error.
        se = sv[:,i]' - sy;
        # Compute the new sx matrix
        sxn = A*sx + se;
        sx = [ sxn; sx[1:order-1,:]];
        # Keep track of some statistics.
        sum_se2 = sum_se2 + sum((se .- sum(se)/length(se)).^2);
        max_sx = maximum([max_sx abs.(sxn)]);
        max_sy = maximum([max_sy abs.(sy)]);
    end
    sigma_se = sqrt(sum_se2/(M*N));

    return (sv, sx, sigma_se, max_sx, max_sy)
end