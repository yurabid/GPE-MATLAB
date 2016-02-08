function [resp,resm] = detect_core(slice,rx,ry)
[n, m] = size(slice);
dx = rx(1,2)-rx(1,1);
dy = ry(2,1)-ry(1,1);
an = angle(slice);
phi = abs(del2(abs(slice).^2,dx,dy));
%idx = [-1 -m 1 m -1];
%idx = [-m-1 -m+1 m+1 m-1 -m-1];
%idx = [0 1 m+1 m 0];
idx = [0 0; 1 0; 1 1; 0 1];
nshift = size(idx,1);
ans = zeros(n,m,nshift);
dif = zeros(n,m,nshift-1);
for i=1:nshift
	ans(:,:,i) = circshift(an,idx(i,:));
end
dif = ans - circshift(ans,[0 0 1]);

res1 = (sum(dif>pi,3)-sum(dif<-pi,3)).*(phi>0.00001);

%an1 = circshift(an,[1,0]);
%an2 = circshift(an,[1,1]);
%an3 = circshift(an,[0,1]);

%diff1 = an1 - an;
%diff2 = an2 - an1;
%diff3 = an3 - an2;
%diff4 = an - an3;

%res = ((diff1>pi) + (diff2>pi) + (diff3>pi) + (diff4>pi) - (diff1<-pi) - (diff2<-pi) - (diff3<-pi) - (diff4<-pi))>0;

res = res1>0;
resp = [rx(res)-dx/2 ry(res)-dy/2];
res = res1<0;
resm = [rx(res)-dx/2 ry(res)-dy/2];
