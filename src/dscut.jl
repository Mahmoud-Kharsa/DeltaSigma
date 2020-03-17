function dscut(p1,y1,p2,y2)
    k1 = y2./(y2-y1)
    k2 = 1-k1
    n = size(p1,1)
    x1=ones(Float64,n,1)
    p = (k1*x1).*p1 + (k2*x1).*p2
    return p
end