function partitionABCD(ABCD, m)
# Partition ABCD into A, B, C, D for an m-input state-space system.

    n = size(ABCD,2)-m; 
    r = size(ABCD,1)-n;

    A = ABCD[1:n, 1:n];
    B = ABCD[1:n, n+1:n+m];
    C = ABCD[n+1:n+r, 1:n];
    D = ABCD[n+1:n+r, n+1:n+m];
	println(A)
	println(B)
	println(C)
	println(D)
	return A, B, C, D
	end