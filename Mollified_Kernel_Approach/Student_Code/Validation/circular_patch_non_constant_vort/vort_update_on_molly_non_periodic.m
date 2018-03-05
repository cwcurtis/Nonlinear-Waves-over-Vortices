function u = vort_update_on_molly_non_periodic(mu,gam,ep,F,u,gvals,Nvorts,dt)

k1 = dt*force_terms_on_molly_non_periodic(mu,gam,ep,F,u,gvals,Nvorts);
uh1 = u+k1/2;

k2 = dt*force_terms_on_molly_non_periodic(mu,gam,ep,F,uh1,gvals,Nvorts);
uh2 = u + k2/2;

k3 = dt*force_terms_on_molly_non_periodic(mu,gam,ep,F,uh2,gvals,Nvorts);
uh3 = u+k3;

k4 = dt*force_terms_on_molly_non_periodic(mu,gam,ep,F,uh3,gvals,Nvorts);

u = u + (k1 + 2*(k2+k3) + k4)/6;