%% Bounds with the optimal a parameters for 1st-kind Chebyshev polynomials

% Location of the optimal value
k = 15;
kval = 1:k;
aopt = gen_aopt(k);

figure(1)
subplot(1,2,1);
semilogy(kval,aopt,'rx',...
    kval,(log(kval)./kval).^2,'k--', ...
    kval,(log(kval)./(3*kval)).^2,'k--',"LineWidth",2)
vline(3,'k--');
legend('a_k^*','log(k)^2/k^2')
xticks([0,3,5,10,15])
xlabel("k")

% Entity of the bound
lambdak = zeros(k, 1);
for n=1:k
    lambdak(n) = abs(1/(2*eval_chebcoef1(aopt(k), n)));
end
figure(1)
subplot(1,2,2)
semilogy(kval,lambdak,'rx',kval,1.03*log(kval)./(2*kval.^2),'k--', ...
    kval,log(kval)./(6*kval.^2),'k--',"LineWidth",2)
vline(3,'k--');
legend('\Lambda_k','1.03 \cdot log(k)/{2k^2}')
xticks([0,3,5,10,15])
xlabel("k")