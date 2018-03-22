function nl = direct_solver(rval,xpos,zpos,gvals)
ov2p = 1/(2*pi);

dx = bsxfun(@minus,xpos,xpos');
dz = bsxfun(@minus,zpos,zpos');
dzp = bsxfun(@plus,zpos,zpos');
dx2 = dx.^2;

diff = dx2 + dz.^2 + eye(size(dx));
er = exp(-diff/(2*rval^2));
scv = (1-er).*(1+2*er)./diff;

diffp = dx2 + dzp.^2;
erp = exp(-diffp/(2*rval^2));
scvp = (1-erp).*(1+2*erp)./diffp;

Kx = -dz.*scv + dzp.*scvp;
Kz = dx.*(scv-scvp);

xdot = ov2p*Kx*gvals;
zdot = ov2p*Kz*gvals;

nl = [xdot zdot];