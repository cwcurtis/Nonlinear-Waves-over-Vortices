function u = vort_update_on_molly_fourier(gam,mu,ep,F,u,gvals,Nvorts,Ntrunc,dt)

k1 = dt*force_terms_on_molly_fourier(gam,mu,ep,F,u,gvals,Nvorts,Ntrunc);
uh1 = u+k1/2;

k2 = dt*force_terms_on_molly_fourier(gam,mu,ep,F,uh1,gvals,Nvorts,Ntrunc);
uh2 = u + k2/2;

k3 = dt*force_terms_on_molly_fourier(gam,mu,ep,F,uh2,gvals,Nvorts,Ntrunc);
uh3 = u+k3;

k4 = dt*force_terms_on_molly_fourier(gam,mu,ep,F,uh3,gvals,Nvorts,Ntrunc);

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