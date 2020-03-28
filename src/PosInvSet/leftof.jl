function leftof(p,a,b)
    b = b-a
    a=a'
    x1=ones(Float64,1,size(p,1))
    x2=a*x1
    x2=x2'
    p = p - x2
    y= b[1]*p[:,2] .> b[2]*p[:,1]
    return y
end