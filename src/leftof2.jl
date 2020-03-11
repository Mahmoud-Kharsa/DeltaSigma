function leftof2(p,a,b)
    #function y=leftof(p,a,b)
    #Return 1 if the point p is to the left of the line ab.
    #For a n x 2 list of points p, return a vector of results.
    #
    # Author: Richard Schreier, Oregon State University. <schreier@ece.orst.edu>
    
    #Translate to the origin and do the check
    x1=ones(Float64,1,size(p,1))
    x2=a*x1
    x2=x2'
    p = p - x2
    b = b-a
    y= b[1]*p[:,2] .> b[2]*p[:,1]
    return y
end