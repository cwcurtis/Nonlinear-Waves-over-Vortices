function nl = force_terms_multipole(mu,gam,rval,u,gvals,Nvorts)

xpos = u(1:Nvorts);
zpos = u(Nvorts+1:2*Nvorts);
pval = 15;

Kvec = multi_pole_kernel(xpos,zpos,gvals,gam,rval,pval);
Kx = Kvec(:,1);
Kz = Kvec(:,2);

xdot = mu/(2*pi)*Kx;
zdot = mu/(2*pi*gam)*Kz;

nl = [xdot;zdot];
