function a_opt = gen_aopt(k)
%GEN_AOPT Generates the optimal for the definition of the 1st-kind
%Chebyshev polynomial accelarted smoothers.
a_opt = zeros(k,1);

for l=1:k
    y = fzero(@(x) aopt_fun(x, l),.5);
    a_opt(l) = y^2;
end



end

function y = aopt_fun(x, n)
% Evaluate the function (1-x^2)^(2n)[(8n-2x)(1-x)^(2n) - (8n+2x)(1+x)^(2n)] + x[(1+x)^(2n) + (1-x)^(2n)][(1+x)^(4n) + (1-x)^(4n)]
%(hopefully) in a stable manner
p_2n_sum = zeros(2*n - 1, 1);
p_2n_diff = zeros(2*n+1, 1);
p_4n = zeros(4*n - 1, 1);
for j = 0:n
    p_2n_sum(2*j+1) = gamma(2*n+1)/(gamma(2*j+1)*gamma(2*n-2*j+1));
    if j < n
        p_2n_diff(2*j+2) = gamma(2*n+1)/(gamma(2*j+2)*gamma(2*n-2*j));
    end
end
for j = 0:2*n
    p_4n(2*j+1) = gamma(4*n+1)/(gamma(2*j+1)*gamma(4*n-2*j+1));
end
qt1 = 2*polyval(p_2n_sum(end:-1:1), x);
qt2 = -2*polyval(p_2n_diff(end:-1:1), x);
qt3 = 2*polyval(p_4n(end:-1:1), x);

y = (1 + x).^(2*n) .* (1 - x).^(2*n) .* (8 * n * qt2 - 2* x.* qt1) + x .* qt1 .* qt3;

end
