function Kmatx_f = kernel_tester(Nvorts,gam,ep,Ntrunc)

xpos = linspace(-.9,.9,Nvorts)';
zpos = linspace(.01,.99,Nvorts)';

Kmatx_f = zeros(Nvorts);
Kmatz_f = zeros(Nvorts);

parfor jj = 1:Nvorts    
    dx = xpos(jj) - xpos;
    dzm = zpos(jj) - zpos;
    
    % compute the mollified part of the kernel
    Kmxm = zeros(Nvorts,1);
    Kmzm = zeros(Nvorts,1);
    
    %spx = sin(pi*dx);
    %cpx = cos(pi*dx);
    
    %facm = 4*(cosh(gam*pi*dzm) - cpx);
    
    %kerxm = -sinh(gam*pi*dzm)./facm;
    %kerzm = spx./facm;
    %kerxm(jj) = 0;
    %kerzm(jj) = 0;
        
    for ll=1:Nvorts                        
        [Kmxm(ll),Kmzm(ll)] = kernel_mol(dx(ll),dzm(ll),gam,ep,Ntrunc);
    end
    
    %Kmatx_f(jj,:) = kerxm + Kmxm;
    %Kmatz_f(jj,:) = kerzm + Kmzm;     
    
    Kmatx_f(jj,:) = Kmxm;
    Kmatz_f(jj,:) = Kmzm;     
end