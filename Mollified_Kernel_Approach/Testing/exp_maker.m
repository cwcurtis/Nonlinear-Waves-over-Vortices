function [kxexp,kzexp] = exp_maker(xloc,zloc,xfar,zfar,gam)

pg = pi*gam;

tl = tanh(pg*zfar);
sl = sech(pg*zfar);
ctl = cos(pi*xfar);
stl = sin(pi*xfar);
c2tl = cos(2*pi*xfar);
s2tl = sin(2*pi*xfar);

tsl = tl.*sl;

t = tanh(pg*zloc);
s = sech(pg*zloc);
ct = cos(pi*xloc);
st = sin(pi*xloc);
c2t = cos(2*pi*xloc);
s2t = sin(2*pi*xloc);

ts = t*s;
t2s = ts*t;
t2 = t^2;

kxexp = .5*(s^2*tl + t^2*tl.^3 + .5*s^2*sl.^2.*tl + ct*ctl.*tl/2 - 2*s*t^2*(ct*sl.*tl.*ctl + st*sl.*tl.*ctl) + st*stl.*tl/2 + .5*s^2*(c2t*sl.^2.*tl.*c2tl + s2t*sl.^2.*tl.*s2tl));        
kzexp = .5*( ts*(st*(ctl.*tsl)-ct*(stl.*tsl)) + s*ts*(s2t*(tsl.*c2tl) - c2t*(tsl.*s2tl)) );

