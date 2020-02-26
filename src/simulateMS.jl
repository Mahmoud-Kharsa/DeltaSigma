using ControlSystems
using MAT
using Polynomials

function simulateMS(v, M, mtf, d, dw, sx0)

    dw = []
    sx0 = []
    M = 16
    d = 0 

    # #Case 1
    # file = matopen("v1.mat")
    # v = read(file, "v")
    # close(file)
    # mtf = zpk([1], [0], 1, 1)

    #Case 2
    # file = matopen("v2.mat")
    # v = read(file, "v")
    # close(file)
    # mtf = zpk([ 1, 1 ], [ 0.3, 0.3 ], 1, 1)
    
    #Case 3
    # file = matopen("v3.mat")
    # v = read(file, "v")
    # close(file)
    # mtf = zpk([1], [0], 1, 1)

    #Case 4
    # file = matopen("v4.mat")
    # v = read(file, "v")
    # close(file)
    # mtf = zpk([1], [0], 1, 1)
    # d = 0.5

    #Case 5
    # file = matopen("v5.mat")
    # v = read(file, "v")
    # close(file)
    # mtf = synthesizeNTF(2,25,1,2)

    #Case 6
    file = matopen("v6.mat")
    v = read(file, "v")
    close(file)
    mtf = synthesizeNTF(4,25*0.9,1,1.3)



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