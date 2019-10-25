function x = bar_pos(t,xm,tm)

if(tm==0)
    x=xm;
else
    delta = tm/7;
    a = xm/(2*delta*(tm-delta));
    x = a*t.^2.*(t<=delta) + (2*a*delta*t - a*delta^2).*(t>delta&t<(tm-delta)) + (xm-a*(tm-t).^2).*(t>=(tm-delta)&t<tm) + xm.*(t>=tm);
%     a=2*xm/tm^2;
%     x = a*t.^2/2.*(t<tm) + xm.*(t>=tm);
end

end
