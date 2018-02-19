function nl = force_terms_on_molly_fourier(gam,mu,ep,F,u,gvals,Nvorts,Ntrunc)

xpos = u(1:Nvorts);
zpos = u(Nvorts+1:2*Nvorts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ntrg = Nvorts*(Nvorts-1)/2;
Nproc = feature('numcores');
Nstep = floor(Ntrg/Nproc);
Nrem = mod(Ntrg,Nproc);
Kmmx = zeros(Nvorts);
Kz = zeros(Nvorts);
Kmvlsx = zeros(Nstep,Nproc);
Kvlsz = zeros(Nstep,Nproc);

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
%ticBytes(gcp);
dxmm = xp1m - xp2m;
dzmm = zp1m - zp2m;

for kk = 1:Nproc    
    dx = dxmm(:,kk);
    dzm = dzmm(:,kk);
    
    spx = sin(pi*dx);
    cpx = cos(pi*dx);
    
    facm = 4*(cosh(gam*pi*dzm) - cpx);
    
    kerxm = -sinh(gam*pi*dzm)./facm;
    kerzm = spx./facm;
    
    [Kmxm,Kmzm] = kernel_mol(dx,dzm,gam,ep,Ntrunc);
    
    Kmvlsx(:,kk) = kerxm + Kmxm;
    Kvlsz(:,kk) = kerzm + Kmzm;
end
%tocBytes(gcp);
Kmvlsx = Kmvlsx(:);
Kvlsz = Kvlsz(:);

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
    Kvlsz = [Kvlsz;(kerzm+Kmzm)];
end

Nprior = 0;
for ll=1:Nvorts-1
   svec = (Nprior + 1):(Nprior + Nvorts-ll); 
   Kmmx = Kmmx + diag(Kmvlsx(svec),ll);
   Kz = Kz + diag(Kvlsz(svec),ll);
   Nprior = Nprior + Nvorts-ll;   
end
    
Kx = Kmmx - Kmmx';
Kz = Kz - Kz';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xdot = mu*F*Kx*gvals;
zdot = mu*F/gam*Kz*gvals;

nl = [xdot;zdot];