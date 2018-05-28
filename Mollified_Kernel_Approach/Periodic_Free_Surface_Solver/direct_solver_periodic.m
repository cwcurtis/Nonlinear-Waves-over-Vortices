function nl = direct_solver_periodic(rval,xpos,zpos,gvals,Mx)
p2M = pi/(2*Mx);

dx = bsxfun(@minus,xpos,xpos');
dz = bsxfun(@minus,zpos,zpos');
dzp = bsxfun(@plus,zpos,zpos');

dzz = dx + 1i*dz + eye(size(dx));
dzzp = dx + 1i*dzp;

diff = abs(dzz).^2;
diffp = abs(dzzp).^2;

er = exp(-diff/(2*rval^2));
scv = (1-2*er).*er./(p2M*dzz);
erp = exp(-diffp/(2*rval^2));
scvp = (1-2*erp).*erp./(p2M*dzzp);

pker = cot(p2M*dzz) + scv;
pker = pker - diag(diag(pker));
pkerp = cot(p2M*dzzp) + scvp;
ptot = pker - pkerp;
Kvals = ptot*gvals;

nl = [imag(Kvals) real(Kvals)];