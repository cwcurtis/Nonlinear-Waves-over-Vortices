function nl = force_terms_on_molly_non_periodic(mu,gam,ep,u,gvals,Nvorts)

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
    dzm = gam*dzmm(:,kk);
    
    rvals = dx.^2 + dzm.^2;
    er = exp(-rvals/(2*ep^2));
    scv = (1-er).*(1+2*er)./rvals;
    kerxm = -dzm.*scv;
    kerzm = dx.*scv;
    
    Kmvlsx(:,kk) = kerxm;
    Kvlsz(:,kk) = kerzm;
end
%tocBytes(gcp);
Kmvlsx = Kmvlsx(:);
Kvlsz = Kvlsz(:);

if Nrem>0   
    dx = xpos(inds(Nproc*Nstep+1:Ntrg,1)) - xpos(inds(Nproc*Nstep+1:Ntrg,2));
    dzm = gam*(zpos(inds(Nproc*Nstep+1:Ntrg,1)) - zpos(inds(Nproc*Nstep+1:Ntrg,2)));
    
    rvals = dx.^2 + dzm.^2;
    er = exp(-rvals/(2*ep^2));
    scv = (1-er).*(1+2*er)./rvals;
    kerxm = -dzm.*scv;
    kerzm = dx.*scv;
    
    Kmvlsx = [Kmvlsx;kerxm];
    Kvlsz = [Kvlsz;kerzm];
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

xdot = mu/(2*pi)*Kx*gvals;
zdot = mu/(2*pi*gam)*Kz*gvals;

nl = [xdot;zdot];