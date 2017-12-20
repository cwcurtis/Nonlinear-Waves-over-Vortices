function Kloc = near_neighbor_comp(xloc,zloc,xfar,zfar,gnear,gfar,Ntrunc,gam,ep)

nparts = length(xloc);
nnghbs = length(xfar);
Kloc = zeros(nparts,2);

for jj=1:nparts
    
    dx = xloc(jj) - xfar;
    dzm = zloc(jj) - zfar;
    dzp = zloc(jj) + zfar;
    spx = sin(pi*dx);
    cpx = cos(pi*dx);
    
    facm = 4*(cosh(gam*pi*dzm) - cpx);
    facp = 4*(cosh(gam*pi*dzp) - cpx);
    
    kerxm = -sinh(gam*pi*dzm)./facm;
    kerzm = spx./facm;
    
    kerxp = -sinh(gam*pi*dzp)./facp;
    kerzp = spx./facp;
    
    % compute the mollified part of the kernel
    Kmx = zeros(nnghbs,1);
    Kmz = zeros(nnghbs,1);

    for ll=1:nnghbs                        
        [Kmx(ll),Kmz(ll)] = kernel_mol(dx(ll),dzm(ll),gam,ep,Ntrunc);                
    end
    
    Kloc(jj,1) = sum(gfar.*(kerxm+Kmx-kerxp));
    Kloc(jj,2) = sum(gfar.*(kerzm+Kmz-kerzp));     
    
end    

for jj=1:nparts
    
    dx = xloc(jj) - xloc;
    dzm = zloc(jj) - zloc;
    dzp = zloc(jj) + zloc;
    spx = sin(pi*dx);
    cpx = cos(pi*dx);
    
    facm = 4*(cosh(gam*pi*dzm) - cpx);
    facp = 4*(cosh(gam*pi*dzp) - cpx);
    
    kerxm = -sinh(gam*pi*dzm)./facm;
    kerzm = spx./facm;
    
    kerxp = -sinh(gam*pi*dzp)./facp;
    kerzp = spx./facp;
    
    kerxm(jj) = 0;
    kerzm(jj) = 0;
    
    % compute the mollified part of the kernel
    Kmx = zeros(nparts,1);
    Kmz = zeros(nparts,1);

    for ll=1:nparts                        
        [Kmx(ll),Kmz(ll)] = kernel_mol(dx(ll),dzm(ll),gam,ep,Ntrunc);                
    end
    
    Kloc(jj,1) = Kloc(jj,1) + sum(gnear.*(kerxm+Kmx-kerxp));
    Kloc(jj,2) = Kloc(jj,2) + sum(gnear.*(kerzm+Kmz-kerzp));     
    
end    