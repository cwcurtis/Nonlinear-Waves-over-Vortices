function nl = force_terms_on_molly_non_periodic(u,omega,gval,Nvorts)

xpos = u(1:Nvorts);
zpos = u(Nvorts+1:2*Nvorts);

dx = bsxfun(@minus,xpos,xpos');
dz = bsxfun(@minus,zpos,zpos');
diff = dx.^2 + dz.^2 + eye(size(dx));
Kx = dz./diff;
Kz = -dx./diff;
xdot = gval/(2*pi)*sum(Kx,2);
zdot = gval/(2*pi)*sum(Kz,2);
nl = [-omega*zpos + xdot;omega*xpos + zdot];
%nl = [xdot;zdot];