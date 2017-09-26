function Kloc = near_neighbor_comp(xloc,zloc,xfar,zfar,gnear,gfar,Ntrunc,gam,ep)

nparts = length(xloc);
nnghbs = length(xfar);
Kloc = zeros(nparts,2);

for jj=1:nparts
    xpjj = xloc(jj);
    zpjj = zloc(jj);

    % Compute out-of-cell nearest neighbors
    for ll=1:nnghbs  
        dx = xpjj - xfar(ll);
        dzm = zpjj - zfar(ll);
        dzp = zpjj + zfar(ll);
        
        [msumxm,msumzm] = kernel_mol(dx,dzm,gam,ep,Ntrunc);
        
        fac = cosh(gam*pi*dzp) - cos(pi*dx);
        msumxp = -sinh(gam*pi*dzp)/fac;
        msumzp = sin(pi*dx)/fac;
        
        Kxm = msumxm;
        Kxp = msumxp;
        Kz = msumzm-msumzp;
        Kloc(jj,1) = Kloc(jj,1) + gfar(ll)*(Kxm-Kxp);
        Kloc(jj,2) = Kloc(jj,2) + gfar(ll)*Kz;
    end    
    
end    

% Compute in-cell nearest neighbors
Kxm = zeros(nparts);
Kxp = zeros(nparts);
Kz = zeros(nparts);

for jj=2:nparts
    xpjj = xloc(jj);
    zpjj = zloc(jj);
    for ll=1:jj-1        
        dx = xpjj - xloc(ll);
        dzm = zpjj - zloc(ll);
        dzp = zpjj + zloc(ll);
        
        [msumxm,msumzm] = kernel_mol(dx,dzm,gam,ep,Ntrunc);
        
        fac = cosh(gam*pi*dzp) - cos(pi*dx);
        msumxp = -sinh(gam*pi*dzp)/fac;
        msumzp = sin(pi*dx)/fac;
        
        Kxm(ll,jj) = msumxm;
        Kxp(ll,jj) = msumxp;
        Kz(ll,jj) = msumzm-msumzp;
    end
end

Kxm = Kxm - Kxm';
Kz = Kz - Kz';
Kxp = Kxp + Kxp';

   
for ll=1:nparts
    Kloc(ll,1) = 1/tanh(gam*pi*zloc(ll)) + sum(gnear.*(Kxm(:,ll)-Kxp(:,ll)));
    Kloc(ll,2) = sum(gnear.*Kz(:,ll));
end