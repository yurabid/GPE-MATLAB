function rho=dens_states(e,mu,om)

ap = 2*(e+mu)/om;
am = 2*(mu-e)/om;

x=sqrt(2*mu/om);
um=sqrt(2*e/om);%sqrt(am+x.^2);
Im1 = um.^3.*x/4 + am.*um.*x/8 - am.^2/8.*log(x+um);

x=sqrt(max(0,am));
um=sqrt(-am+x.^2);
Im2 = um.^3.*x/4 + am.*um.*x/8 - am.^2/8.*log(x+um);

x=sqrt(ap);
up=0;%sqrt(ap-x.^2);
Ip1 = -up.^3.*x/4 + ap.*up.*x/8 + ap.^2/8.*pi/2;

x=sqrt(2*mu/om);
up=sqrt(2*e/om);
Ip2 = -up.^3.*x/4 + ap.*up.*x/8 + ap.^2/8.*real(asin(x./sqrt(ap)));

rho = 2/pi/om*real(Im1-Im2 + Ip1-Ip2);
