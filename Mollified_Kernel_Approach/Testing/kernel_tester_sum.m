function [Kmatx_s] = kernel_tester_sum(Nvorts,gam,ep,Mterms)

xpos = linspace(-.9,.9,Nvorts);
zpos = linspace(.01,.99,Nvorts);

%efun = @(x,z,m) (1-exp(-((x-2.*m).^2+gam.^2*z.^2)/ep^2))./(2*pi*( (x-2.*m).^2 + gam.^2.*z.^2 ));
efun = @(x,z,m) (-exp(-((x-2.*m).^2+gam.^2*z.^2)/ep^2))./(2*pi*( (x-2.*m).^2 + gam.^2.*z.^2 ));
Mvec = -Mterms:Mterms;

Kmatx_s = zeros(Nvorts);
Kmatz_s = zeros(Nvorts);

for jj = 1:Nvorts    
    dx = xpos(jj) - xpos;
    dzm = zpos(jj) - zpos;
    
     % compute the mollified part of the kernel
    Kmxs = zeros(Nvorts,1);
    Kmzs = zeros(Nvorts,1);
    
    %spx = sin(pi*dx);
    %cpx = cos(pi*dx);
    
    %facm = 4*(cosh(gam*pi*dzm) - cpx);
    
    %kerxm = -sinh(gam*pi*dzm)./facm;
    %kerzm = spx./facm;
    %kerxm(jj) = 0;
    %kerzm(jj) = 0;
        
    for ll=1:Nvorts           
        if dx(ll)~=0 || dzm(ll) ~=0
            evec = efun(dx(ll),dzm(ll),Mvec);
            Kmxs(ll) = -gam*dzm(ll)*sum(evec);
            Kmzs(ll) = sum((dx(ll)-2*Mvec).*evec);            
            plot(Mvec,evec)
            pause
        end         
    end
    
    Kmatx_s(jj,:) = Kmxs;
    Kmatz_s(jj,:) = Kmzs;
end






