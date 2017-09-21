function msumz = kernel_vert(x,z,gam,ep,Ntrunc)

    Mvals = 1:Ntrunc;

    ghat=@(z,m) 1i/4*(exp(gam*z*pi*m).*(-E1(z,m,gam,ep)/2)+exp(-gam*z*pi*m).*(E2(z,m,gam,ep)/2-1));
    
    if x == 0
    
        msumz = 0;
        
    else
        
        msumz = -2*sum( sin(pi*Mvals*(x)).*imag( ghat(z,Mvals)) );
        
    end
    