function nl = force_terms_on_molly_non_periodic(mu,rval,u,gvals,Nvorts)

xpos = u(1:Nvorts);
zpos = u(Nvorts+1:2*Nvorts);

dx = bsxfun(@minus,xpos,xpos');
dz = bsxfun(@minus,zpos,zpos');
diff = dx.^2 + dz.^2 + eye(size(dx));
er = exp(-diff/(2*rval^2));
scv = (1-er).*(1+2*er)./diff;
Kx = (-dz.*scv);
Kz = (dx.*scv);
xdot = mu/(2*pi)*Kx*gvals;
zdot = mu/(2*pi)*Kz*gvals;
nl = [xdot;zdot];