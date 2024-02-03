function u = vort_update_on_molly_non_periodic(u,omega,gval,Nvorts,dt)

k1 = dt*force_terms_on_molly_non_periodic(u,omega,gval,Nvorts);
uh1 = u+k1/2;

k2 = dt*force_terms_on_molly_non_periodic(uh1,omega,gval,Nvorts);
uh2 = u + k2/2;

k3 = dt*force_terms_on_molly_non_periodic(uh2,omega,gval,Nvorts);
uh3 = u+k3;

k4 = dt*force_terms_on_molly_non_periodic(uh3,omega,gval,Nvorts);

u = u + (k1 + 2*(k2+k3) + k4)/6;