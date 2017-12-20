function [msumx,msumz] = kernel_mol(x,z,gam,ep,Ntrunc)

    if x ~= 0 || z ~=0
        alpha = gam*abs(z)/ep;
        beta = gam*pi*abs(z);
        sz = sign(z);
        fac1 = @(m) 1/4*exp(-beta*m).*erfc(alpha - ep*pi*m/2);
        fac2 = @(m) 1/4*exp(beta*m).*erfc(alpha + ep*pi*m/2);
        %fac3 = @(m) 1/2*exp(-beta*m);
        
        Mvals = 1:Ntrunc;
        f1 = fac1(Mvals);
        f2 = fac2(Mvals);
        
        %umol = fac3(Mvals);
        %fhat = -umol + f1 + f2;
        %ghat = umol - (f1 - f2);    
        %msumx = sz*(-1/4 + fac1(0)/2 + sum( cos(pi*Mvals*x).*fhat ) );
        %msumz = sum( sin(pi*Mvals*x).*ghat );
        
        fhat = f1 + f2;
        ghat = -(f1 - f2);    
        msumx = sz*( fac1(0)/2 + sum( cos(pi*Mvals*x).*fhat ) );
        msumz = sum( sin(pi*Mvals*x).*ghat );
    else
        msumx = 0;
        msumz = 0;
    end
    