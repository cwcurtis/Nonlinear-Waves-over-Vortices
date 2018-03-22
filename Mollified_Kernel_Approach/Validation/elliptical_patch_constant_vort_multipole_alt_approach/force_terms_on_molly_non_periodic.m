function nl = force_terms_on_molly_non_periodic(mu,gam,rval,u,gvals,Nvorts)
m2p = mu/(2*pi);

xpos = u(1:Nvorts);
zpos = u(Nvorts+1:2*Nvorts);

dx = bsxfun(@minus,xpos,xpos');
dz = gam*bsxfun(@minus,zpos,zpos');
dzp = gam*bsxfun(@plus,zpos,zpos');
dx2 = dx.^2;

diff = dx2 + dz.^2 + eye(size(dx));
er = exp(-diff/(2*rval^2));
scv = (1-er).*(1+2*er)./diff;

diffp = dx2 + dzp.^2;
erp = exp(-diffp/(2*rval^2));
scvp = (1-erp).*(1+2*erp)./diffp;

Kx = -dz.*scv + dzp.*scvp;
Kz = dx.*(scv-scvp);

xdot = m2p*Kx*gvals;
zdot = m2p/gam*Kz*gvals;

nl = [xdot;zdot];