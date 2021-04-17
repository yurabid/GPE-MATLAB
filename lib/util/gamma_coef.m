function k = gamma_coef(mu,ecut,temp)
% Computes gamma coefficient for SPGPE
k=0;
beta = 1/temp;
for i=1:100
	k = k + exp(((mu-2*ecut)*i+mu)*beta) * lerch(exp((mu-ecut)*beta),1,i)^2;
end