include("dscut.jl")
function dssplit2d(u,ABCD,p)
    n = size(ABCD,1)-1
    C = ABCD[n+1,1:n]'
    D1= ABCD[n+1, n+1]
    N = size(p,2)
    if length(u)==1
        D1u = D1*u
        y = C*p + D1u*ones(Float64,1,N)
    else
        if D1 != 0
            println("D1 must be zero when u is a range")
            return
        else
            y = C*p
        end
    end
    sign1=sign(sign(y[1])+.5)
    i=findall(y->sign(sign(y)+.5)!=sign1,y)
    i1 = i[1][2]
    pa = dscut(p[:,i1-1],y[i1-1],p[:,i1],y[i1])
    l=length(i)
    i2 = i[l][2]
    pb = dscut(p[:,i2],y[i2], p[:,i2+1],y[i2+1])
    answer=ones(Int64,length(i),1)
    for il=1:length(i)
        answer[il]=i[il][2]
    end
    if sign1 > 0
        take=p[:,answer[1]:answer[length(answer)]]
        part1=hcat(pa,take)
        part2=hcat(pb,pa)
        pminus=hcat(part1,part2)
        take1=p[:,1:i1-1]
        take2=p[:,i2+1:N]
        pplus = hcat(take1,pa,pb,take2)
    else
        take=p[:,answer[1]:answer[length(answer)]]
        part1=hcat(pa,take)
        part2=hcat(pb,pa)
        pplus=hcat(part1,part2)
        take1=p[:,1:i1-1]
        take2=p[:,i2+1:N]
        pminus = hcat(take1,pa,pb,take2)
    end
    ne = size(pplus,2)
    eplus = vcat(2:ne,1)
    eplus = hcat(1:ne,eplus)
    eplus=eplus'
    ne = size(pminus,2)
    eminus = vcat(2:ne,1)
    eminus = hcat(1:ne,eminus)
    eminus=eminus'
    return pplus,eplus,pminus,eminus
end