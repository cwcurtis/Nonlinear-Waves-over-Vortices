function [msumx,msumz] = kernel_mol(x,z,gam,ep,Ntrunc)

    Mvals = ones(length(x),1)*(1:Ntrunc);
    
    xvals = x*ones(1,Ntrunc);
    zvals = z*ones(1,Ntrunc);
    
    alpha = gam*abs(zvals)/ep;
    beta = gam*pi*abs(zvals);
    
    sz = sign(z);
    emat = exp(-alpha.^2).*exp(-(pi*Mvals*ep/2).^2);
    
    kxzm = -1/(2*sqrt(pi))*alpha.*exp(-alpha.^2);
    kxho = -1/(2*sqrt(pi))*alpha.*emat;
    kzho = -1/(2*sqrt(pi))*(pi*ep*Mvals/2).*emat;
    
    fac1 = 1/4*exp(-beta.*Mvals).*erfc(alpha - ep*pi*Mvals/2);
    fac2 = 1/4*exp(beta.*Mvals).*erfc(alpha + ep*pi*Mvals/2);
    
    fac10 = 1/4*erfc(alpha);
    
    fmat = fac1+fac2;
    gmat = fac1-fac2;
    
    msumx = sz.*sum( (fac10+kxzm)/2 + cos(pi*Mvals.*xvals).*(fmat+kxho),2);
    msumz = -sum( sin(pi*Mvals.*xvals).*(gmat+kzho),2);
    