function u = vort_update_on_molly_non_periodic(mu,gam,ep,u,gvals,Nvorts,dt)

k1 = dt*force_terms_on_molly_non_periodic(mu,gam,ep,u,gvals,Nvorts);
uh1 = u+k1/2;

k2 = dt*force_terms_on_molly_non_periodic(mu,gam,ep,uh1,gvals,Nvorts);
uh2 = u + k2/2;

k3 = dt*force_terms_on_molly_non_periodic(mu,gam,ep,uh2,gvals,Nvorts);
uh3 = u+k3;

k4 = dt*force_terms_on_molly_non_periodic(mu,gam,ep,uh3,gvals,Nvorts);

u = u + (k1 + 2*(k2+k3) + k4)/6;

xpos = u(1:Nvorts);
zpos = u(Nvorts+1:2*Nvorts);

% Relocate vortices that have crossed the boundary

inds = xpos < -1;
xpos(inds) = xpos(inds) + 2;
inds = xpos > 1;
xpos(inds) = xpos(inds) - 2;

u(1:Nvorts) = xpos;
u(Nvorts+1:2*Nvorts) = zpos;