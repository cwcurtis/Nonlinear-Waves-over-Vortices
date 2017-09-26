function Kfar = far_panel_comp(xloc,zloc,xfar,zfar,gfar,gam)

nparts = length(xloc);
Kfar = zeros(nparts,2);
pg = pi*gam;

tl = tanh(pg*zfar);
sl = sech(pg*zfar);
ctl = cos(pi*xfar);
stl = sin(pi*xfar);
ct2l = cos(2*pi*xfar);
st2l = sin(2*pi*xfar);

tsl = gfar.*tl.*sl;
ts2l = tsl.*sl;
t2l = gfar.*tl.^2;
t2sl = t2l.*sl;

t = tanh(pg*zloc);
s = sech(pg*zloc);
c = cosh(pg*zloc);
ct = cos(pi*xloc);
st = sin(pi*xloc);
ct2 = cos(2*pi*xloc);
st2 = sin(2*pi*xloc);

s2 = s.^2;
ts = t.*s;
ts2 = ts.*s;

Kfar(:,1) = s2*sum(gfar.*tl.^3)...
            +(2*s-c).*s2.*(ct*sum(tsl.*ctl)+st*sum(tsl.*stl)) ...
            -3*s2.*st2*sum(ts2l.*st2l)/2-s2*sum(ts2l.*(ct2l-st2l))/2 ...
            -s2.*ct2*sum(ts2l.*ct2l)/2;
        
Kfar(:,2) = ts.*st.*sum(t2l.*ctl)-ts.*ct.*sum(t2l.*stl)...
            +ts2.*st2*sum(t2sl.*ct2l)-ts2.*ct2*sum(t2sl.*st2l);

end

