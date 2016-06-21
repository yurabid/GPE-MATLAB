term = load('term_1p_new.txt','-ascii');
%s = fitoptions('Method','NonlinearLeastSquares','Robust','on','Lower',[0,0,0,0,0,0,0,0]);
%s = fitoptions('Method','NonlinearLeastSquares','Robust','on','StartPoint',[14.54, 33.96, 14.04, 27.55, -30.28, 31.3, -12.21, 2.278]);
weights=zeros(1,size(term,1));
weights(:)=0.01;
weights(45:75)=1000;
s = fitoptions('Method','NonlinearLeastSquares','Robust','on','StartPoint',[1,1,1,1,1,1,1,1],'Weights',weights(:),'MaxFunEvals',10000,'Lower',[0,0,0,0,0,0,0,0]);
%f = fittype('a*(x+b)^n/(x+c)^m','options',s);
%f = fittype('(1+a*x)/(b+c*x+d*x^2+e*x^3+32/3*a*x^4)');
%f = fittype('-(-1+x*(6*1-a1)+x^2*(6*a1-a2)+x^3*(3/32*a5+21/256*a6+135/2048*a7+159/1024*a8)+x^4*(3/32*a6+21/256*a7+135/2048*a8)+x^5*(3/32*a7+21/256*a8)+x^6*(3/32*a8))/(1*x+a1*x^2+a2*x^3+a3*x^4+a4*x^5+a5*x^6+a6*x^7+a7*x^8+a8*x^9)','options',s);
f = fittype('-(-1+x*(-10/9*1-a1)+x^2*(-10/9*a1-a2)+x^3*(3/32*a5+21/256*a6+135/2048*a7+159/1024*a8)+x^4*(3/32*a6+21/256*a7+135/2048*a8)+x^5*(3/32*a7+21/256*a8)+x^6*(3/32*a8)+32/pi*exp(-2*x-1)*x^9*a8)/(1*x+a1*x^2+a2*x^3+a3*x^4+a4*x^5+a5*x^6+a6*x^7+a7*x^8+a8*x^9)','options',s);
%f = fittype('-(-1+x*(-10/9*1-a1)+x^2*(-10/9*a1-a2)+x^3*(3/32*a5+21/256*a6+135/2048*a7)+x^4*(3/32*a6+21/256*a7)+x^5*(3/32*a7)+32/pi*exp(-2*x-1)*x^8*a7)/(1*x+a1*x^2+a2*x^3+a3*x^4+a4*x^5+a5*x^6+a6*x^7+a7*x^8)','options',s);
%f = fittype('-(-1+x*(6*1-a1)+x^2*(6*a1-a2)+x^3*(3/32*a5+21/256*a6+135/2048*a7)+x^4*(3/32*a6+21/256*a7)+x^5*(3/32*a7))/(1*x+a1*x^2+a2*x^3+a3*x^4+a4*x^5+a5*x^6+a6*x^7+a7*x^8)','options',s);
%f = fittype('-(-1+x*(6*1-a1)+x^2*(6*a1-a2)+x^3*(3/32*a5+21/256*a6)+x^4*(3/32*a6))/(1*x+a1*x^2+a2*x^3+a3*x^4+a4*x^5+a5*x^6+a6*x^7)','options',s);
%f = fittype(@(a,b,c,m,n,x) a.*(x+b).^n./(x+c).^m,'options',s);
%[c2,gof2,output] = fit(term(1:end,1),term(1:end,2),f)
[c2,gof2,output] = fit(term(:,1),term(:,2)+1./term(:,1),f)
%[c2,gof2,output] = fit((1:1000)'*0.01,fun((1:1000)'*0.01)+1./((1:1000)'*0.01),f)
figure
hold all
%plot(term(1:end,1),term(1:end,2),'+')
plot(term(1:end,1),term(1:end,2)+1./term(1:end,1),'+')
plot(c2)
%plot((1:1000)*0.01,-c2((1:1000)*0.01).'+fun((1:1000)*0.01)+1./((1:1000)*0.01));
%plot((1:1000)'*0.03,c2((1:1000)'*0.03)+1./((1:1000)'*0.03))
%%
ind = [ 1 2 3 4 5 6 7 ];
%val = [20.40 39.9 63.8 91.6 122.5];
val = [ 6.01 19.03 37.27 59.03 84.43 116.4 155.71];
%val = [55.8 335.5 916.5 1741.5 2807 4156 5766 7919];
%s = fitoptions('polynomial');
f = fittype('a.*(x+1).^2  + c');
%f = fittype('poly2');
[c2,gof2] = fit(ind.',val.',f)
figure
hold all
plot(ind,val,'+')
plot(c2)

%%
ind = [0 1 2];
val = [156 892 2340];
%s = fitoptions('polynomial');
%f = fittype('a.*x.^2 + a.*x + b');
f = fittype('poly2');
[c2,gof2] = fit(ind.',val.',f)
figure
hold all
plot(ind,val,'+')
plot(c2)