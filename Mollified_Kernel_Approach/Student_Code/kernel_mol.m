function [msumx,msumz] = kernel_mol(x,z,gam,ep,Ntrunc)

    if x ~= 0 || z ~=0
        fac1 = @(m) 1/(4*sqrt(pi))*exp(-gam*z*pi*m).*erfc(gam*z/ep - ep*pi*m/2);
        fac2 = @(m) 1/(4*sqrt(pi))*exp(gam*z*pi*m).*erfc(gam*z/ep + ep*pi*m/2);
    
        Mvals = 0:Ntrunc;
        f1 = fac1(Mvals);
        f2 = fac2(Mvals);
        umolfac = .5*exp(-pi*gam*z*Mvals);
    
        fhat = -umolfac + f1 + f2;
        ghat = umolfac + f1 - f2;
    
        msumx = sum( cos(pi*Mvals*x).*fhat );
        msumz = -sum( sin(pi*Mvals*x).*ghat );
    else
        msumx = 0;
        msumz = 0;
    end
    