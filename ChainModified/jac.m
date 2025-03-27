function C = jac(A,B)
    C = max(A,B) + log(1+exp(-abs(A-B)));
end
	