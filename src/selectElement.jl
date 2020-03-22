function selectElement(v, sy, dw, tri=0)

    M = length(sy);

    if isempty(dw) || all(dw .== 1)
        dw_is_one = true;
        sdw = M;
    else
        dw_is_one = false;
        sdw = sum(dw);
    end
    
    if abs(v) > sdw
        error("|v| is too large.");
    elseif v == sdw
        sv = ones(M,1);
        return sv
    elseif v == -sdw
        sv = -ones(M,1);
        return sv
    end
    
    sv = -ones(M,1);
    possibilities = sortperm(vec(-sy)); 	#Determine the element priority
    
    # Speed up selection in the usual case where all element weights are one.
    # Suggested 1997-03-05 by J. A. Cherry.
    if dw_is_one
        sv[possibilities[1:floor(Int, (M+v)/2)]] .= 1;
        return sv
    end
    
    # Go through sv possibilities one by one, until one which meets the
    # v = sv'*dw constraint is found.
    error("dw not ones has a bug currently");
    
    i = 1;				    # Selection level
    pointer =  ones(1,M);	# Array of pointers to selected elements
    selected = zeros(1,M);	# Selected elements
    dw0 = [0; dw];
    while true
        while pointer[i] > M
            # backtrack
            i = i - 1;
            if i == 0
                break;	#failure!
            end
            pointer[i] = pointer[i] + 1;
            selected[i+1:end] = 0;
        end
        selected[i] = possibilities[pointer[i]];
        dv = 2 * sum(dw0[selected + 1]) - sdw;
        if dv==v
            break;		#success!
        elseif dv<v
            # Proceed to the next level of selection
            pointer[i+1] = pointer[i]+1;
            i = i + 1;
        else
            # Try the next element at the current level.
            pointer[i] = pointer[i]+1;
        end
    end
    selected = selected[selected .!= 0];
    sv[selected] = 1;    

    return sv

end