using Pkg
Pkg.add("Roots")
using Roots
function evalRPoly(roots,x)
    y = ones(Float64,length(x),1)
    for i=1:length(roots)
        y = y.*(x-roots[i])
    end
end
