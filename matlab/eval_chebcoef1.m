function x = eval_chebcoef1(a, n)
%EVAL_CHEBCOEF1 Returns the coefficient of degree 1 of the chebyshev 
% polynomial of degree n over [a, 1]
	x = sqrt(a);
	num = zeros(2*n - 1, 1);
	denum = zeros(2*n+1, 1);
	for j = 0:n
		if j<n
			num(2*j+1) = nchoosek(2*n, 2*j+1);
		end
		denum(2*j+1) = nchoosek(2*n, 2*j);
	end
	x = -n * polyval(num(end:-1:1), x)./polyval(denum(end:-1:1), x);
	
end
