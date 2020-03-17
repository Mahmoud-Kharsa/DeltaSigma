function outconvex2d(x,p)
    n = size(p,2)
    w1=p[2,2:n] - p[2,1:n-1]
    w2=p[1,1:n-1] - p[1,2:n]
    A = hcat(w1,w2)
    B = zeros(Float64,n-1,1)
    for i=1:n-1
        B[i] = A[i,:]'*p[:,i]
    end
    
    out = sum(A*x .> B*ones(Int64,1,size(x,2)),dims=1)
    return out
end