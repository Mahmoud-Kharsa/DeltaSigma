include("leftof.jl")

function hull2d(p)
    p=p'
    x = minimum(p[:,1])
    leftmost = findall(a->a==x,p[:,1])
    x1=minimum(p[leftmost,2])
    y=minimum(p[leftmost,2], dims = 1)
    indices=findall(b->b==y[1],p)
    il=indices[1][1]
    l = [x y]

    x = maximum(p[:,1])
    rightmost = findall(a->a==x,p[:,1])
    x1=maximum(p[rightmost,2])
    y=maximum(p[rightmost,2], dims = 1)
    indices=findall(b->b==y[1],p)
    ir=indices[1][1]
    # ir = rightmost[ir]
    r = [x y]

    if l==r		# Degenerate cases
        v = l
        i = il
        return
    end

    #Split the points into those above and those below the l-r line.
    isAbove = leftof(p,l,r)
    # println("isAbove: ",isAbove)
    ia=findall(a->a==1,isAbove)
    ib=findall(a->a==0,isAbove)
    above = p[ia,:]
    below = p[ib,:]

    #Sort them in terms of increasing first coordinate.
    if any(x->x==1,isAbove)
        isort = sortperm(above[:,1])
        first=l
        second=above[isort,:]
        above = vcat(first,second)
        ia = vcat(il,ia[isort])
    else # no points above
        above = l
        ia = il
    end
    isort = sortperm(below[:,1])
    below = vcat(below[isort,:],r)	# includes the l and r points r must be last.
    ib = vcat(ib[isort],ir)

    #Move along the underside, building the vertex list as we go.
    a = below[1,:]
    a=a'
    nb = size(below,1)
    b = below[2,:]
    b=b'
    v = vcat(a,b)
    # println("Get V: ",v)
    i = vcat(ib[1],ib[2])
    nv = 2
    # println("Before: ",size(v))
    count=0
    countone=0
    countzero=0
    for n=3:nb
        # count=count+1
        p = below[n,:]
        # println("Count: ",count)
        # println("A0: ",a)
        # println("B0: ",b)
        get=leftof(p',a,b)
        # println("Value0: ",get)
        if get==Bool[1]
            countone=countone+1
            cond=true
        end
        if get==Bool[0]
            countzero=countzero+1
            cond=false
        end
        while !cond
            nv = nv-1
            # println("1: ",nv)
            v = v[1:nv,:]
            # println("2: ",v)
            i = i[1:nv]
            # println("3: ",i)
            b = a
            # println("4: ",b)
            if nv>1
                a = v[nv-1,:]
                a=a'
                # println("5: ",a)
            else
                break
            end
            get=leftof(p',a,b)
            if get==Bool[1]
                countone=countone+1
                cond=true
            end
            if get==Bool[0]
                countzero=countzero+1
                cond=false
            end
        end
        v = vcat(v,p')
        i = vcat(i,ib[n])
        nv = nv+1
        a = b
        b = p'
        # println("A: ",a)
        # println("B: ",b)
        # println("")
        # println("")
    end
    # println("Count One: ",countone)
    # println("Count Zero: ",countzero)
    # println("After: ",size(v))

    na = size(above,1)
    for n=na:-1:1
        p = above[n,:]
        get=leftof(p',a,b)
        if get==Bool[1]
            countone=countone+1
            cond=true
        end
        if get==Bool[0]
            countzero=countzero+1
            cond=false
        end
        while ~cond && nv>2
            nv = nv-1
            v = v[1:nv,:]
            i = i[1:nv]
            b = a
            a = v[nv-1,:]
            a=a'
            get=leftof(p',a,b)
            if get==Bool[1]
                countone=countone+1
                cond=true
            end
            if get==Bool[0]
                countzero=countzero+1
                cond=false
            end
        end
        v = vcat(v,p')
        i = vcat(i,ia[n])
        nv = nv+1
        a = b
        b = p'
    end
    return v,i
end
