function [xpos,zpos,r0,lscl] = initializer(omega,gval,Nvorts)

rho = 2*omega/gval;
lscl = 1/sqrt(rho);
disp(lscl)

r0 = sqrt(Nvorts/(pi*rho));
disp(r0)
sqN = ceil(2.*sqrt(Nvorts/pi))+1;
xvals = linspace(-r0, r0, sqN);
zvals = xvals;
xpos = [];
zpos = [];

cfun = @(x,z) x.^2 + z.^2;


for jj=1:sqN
    ainds = cfun(xvals(jj),zvals) <= r0.^2;
    flvds = length(zvals(ainds));
    xpert = .1*lscl*2.*(rand(1,flvds) - .5);
    zpert = .1*lscl*2.*(rand(1,flvds) - .5);
    if flvds>=1
        xpos = [xpos xvals(jj)*ones(1,flvds)+xpert];
        zpos = [zpos zvals(ainds)+zpert];        
    end    
end
