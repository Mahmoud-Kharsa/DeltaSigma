function dsmap(u,ABCD,nlev,x,e,v)
    n = size(ABCD,1)-1
    A = ABCD[1:n, 1:n]
    B = ABCD[1:n, n+1:n+2]
    C = ABCD[n+1, 1:n]
    C=C'
    D1= ABCD[n+1, n+1]
    N = size(x,2)

    if length(u)==2
        u2 = u[2]
        u  = u[1]
        isRange = 1
        println("Supportd only 1x1 u value")
    elseif length(u)==1
        isRange = 0
    else
        println("The dimensions of u are wrong")
        return
    end
    
    if length(v)!=N
        v = v*ones(Int64,1,N)
    else
        println("v is wrong")
    end
    
    # 1) Compute the images of the vertices.
    B1u = B[:,1]*u
    x1=ones(Int64,1,N)
    nx = A*x + B1u*x1 + B[:,2]*v
    
    # 2) For those edges which cross a (or several) splitting hyperplanes,
    #    compute the two images of each split point.
    x2=v[e[1,:]]-v[e[2,:]]
    x2=x2'
    for increment=1:length(x2)
        x2[increment]=abs(x2[increment])
    end
    return nx
end
    