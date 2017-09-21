function msumx = kernel_horiz(x,z,gam,ep,Ntrunc)

    fhat= @(s,m) 1/4*(exp(gam*s*pi*m).*E1(s,m,gam,ep)/2+exp(-gam*s*pi*m).*E2(s,m,gam,ep)/2-exp(-gam*s*pi*abs(m)));

    Mvals = 1:Ntrunc;

    msumx = fhat(z,0) + 2*sum( cos(pi*Mvals*x).*real( fhat(z,Mvals) ) );
    