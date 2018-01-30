function nl = force_terms_on_molly_fourier(gam,mu,ep,F,u,gvals,Nvorts,Ntrunc)

xpos = u(1:Nvorts);
zpos = u(Nvorts+1:2*Nvorts);

xdot = zeros(Nvorts,1);
zdot = zeros(Nvorts,1);

Ntrg = Nvorts*(Nvorts-1)/2;
Nproc = 2;
Nstep = floor(Ntrg/Nproc);
Nrem = mod(Ntrg,Nproc);

Kmvlsx = zeros(Nstep,Nproc);
Kmvlsz = zeros(Nstep,Nproc);
Kmmx = zeros(Nvorts);
Kmmz = zeros(Nvorts);
inds = zeros(Ntrg,2);
Nprior = 0;
for jj=1:Nvorts-1
    jvec = (Nprior + 1):(Nprior + Nvorts-jj); 
    inds(jvec,1) = 1:Nvorts-jj;
    inds(jvec,2) = (jj+1):Nvorts;    
    Nprior = Nprior + Nvorts-jj;
end

xp1m = reshape(xpos(inds(1:Nstep*Nproc,1)),Nstep,Nproc);
xp2m = reshape(xpos(inds(1:Nstep*Nproc,2)),Nstep,Nproc);
zp1m = reshape(zpos(inds(1:Nstep*Nproc,1)),Nstep,Nproc);
zp2m = reshape(zpos(inds(1:Nstep*Nproc,2)),Nstep,Nproc);

parfor kk = 1:Nproc
    
    xp1 = xp1m(:,kk);
    xp2 = xp2m(:,kk);
    
    zp1 = zp1m(:,kk);
    zp2 = zp2m(:,kk);
        
    dx = xp1 - xp2;
    dzm = zp1 - zp2;
    
    spx = sin(pi*dx);
    cpx = cos(pi*dx);
    
    facm = 4*(cosh(gam*pi*dzm) - cpx);
    
    kerxm = -sinh(gam*pi*dzm)./facm;
    kerzm = spx./facm;
         
    [Kmxm,Kmzm] = kernel_mol(dx,dzm,gam,ep,Ntrunc);
    
    Kmvlsx(:,kk) = kerxm + Kmxm;
    Kmvlsz(:,kk) = kerzm + Kmzm;
    
end

Kmvlsx = Kmvlsx(:);
Kmvlsz = Kmvlsz(:);

if Nrem>0   
    dx = xpos(inds(Nproc*Nstep+1:Ntrg,1)) - xpos(inds(Nproc*Nstep+1:Ntrg,2));
    dzm = zpos(inds(Nproc*Nstep+1:Ntrg,1)) - zpos(inds(Nproc*Nstep+1:Ntrg,2));
    
    spx = sin(pi*dx);
    cpx = cos(pi*dx);
    
    facm = 4*(cosh(gam*pi*dzm) - cpx);
    
    kerxm = -sinh(gam*pi*dzm)./facm;
    kerzm = spx./facm;
    
    [Kmxm,Kmzm] = kernel_mol(dx,dzm,gam,ep,Ntrunc);

    Kmvlsx = [Kmvlsx;(kerxm+Kmxm)];
    Kmvlsz = [Kmvlsz;(kerzm+Kmzm)];    
end

Nprior = 0;
for ll=1:Nvorts-1
   svec = (Nprior + 1):(Nprior + Nvorts-ll); 
   Kmmx = Kmmx + diag(Kmvlsx(svec),ll);
   Kmmz = Kmmz + diag(Kmvlsz(svec),ll);
   Nprior = Nprior + Nvorts-ll;
end
    
Kmmx = Kmmx - Kmmx';
Kmmz = Kmmz - Kmmz';

for jj=1:Nvorts      
    xdot(jj) = sum(gvals'.*(Kmmx(jj,:)));
    zdot(jj) = sum(gvals'.*(Kmmz(jj,:)));
    
    xdot(jj) = real(mu*F*xdot(jj));
    zdot(jj) = real(mu*F/gam*zdot(jj));        
end

nl = [xdot;zdot];