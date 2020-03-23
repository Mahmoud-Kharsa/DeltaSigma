using Pkg
Pkg.add("Roots")
using Roots
function evalRPoly(ntf_p,z)
    y1 = ones(Complex{Float64},length(z),1)
    for i=1:length(ntf_p)
        y1 = y1.*(z-ntf_p[i])
    end
    return y1[1]
end
