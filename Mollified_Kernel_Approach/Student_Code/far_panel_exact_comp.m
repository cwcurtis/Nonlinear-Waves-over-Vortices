function Kfar = far_panel_exact_comp(xloc,zloc,xfar,zfar,gfar,Ntrunc,gam,ep)

nparts = length(xloc);
nnghbs = length(xfar);
Kfar = zeros(nparts,2);

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

    %for ll=1:nnghbs                        
    %    [Kmx(ll),Kmz(ll)] = kernel_mol(dx(ll),dzm(ll),gam,ep,Ntrunc);                
    %end
    
    %Kfar(jj,1) = sum(gfar.*(kerxm+Kmx-kerxp));
    %Kfar(jj,2) = sum(gfar.*(kerzm+Kmz-kerzp));     
    Kfar(jj,1) = sum(gfar.*(kerxm-kerxp));
    Kfar(jj,2) = sum(gfar.*(kerzm-kerzp));     
    
end    