function expansion_tester(Nvorts,gam)

xpos = linspace(-.9,.9,Nvorts)';
zpos = linspace(.01,.99,Nvorts)';

for jj = 1:Nvorts    
    dx = xpos(jj) - xpos;
    dzm = zpos(jj) - zpos;
    dzp = zpos(jj) + zpos;
    
    spx = sin(pi*dx);
    cpx = cos(pi*dx);
    
    facm = 4*(cosh(gam*pi*dzm) - cpx);
    facp = 4*(cosh(gam*pi*dzp) - cpx);
    
    kerxm = -sinh(gam*pi*dzm)./facm;
    kerzm = spx./facm;
    
    kerxm(jj) = 0;
    kerzm(jj) = 0;
    
    kerxp = -sinh(gam*pi*dzp)./facp;
    kerzp = spx./facp;   
    
    [kxexp,kzexp] = exp_maker(xpos(jj),zpos(jj),xpos,zpos,gam);
    
    figure(1)
    plot(xpos,log10(abs(kerxm-kerxp-kxexp)))
   
    figure(2)
    plot(xpos,log10(abs(kerzm-kerzp-kzexp)))
    pause
end