function u = vort_update_multipole(mu,gam,rval,u,gvals,Nvorts,dt)

pval = 15;

tree_vals = multi_pole_kernel_build(u(1:Nvorts),u(Nvorts+1:2*Nvorts),gvals,gam,rval,pval);
k1 = dt*force_terms_multipole(mu,gam,rval,u,gvals,Nvorts,pval,tree_vals);
uh1 = u+k1/2;

tree_vals = multi_pole_kernel_build(uh1(1:Nvorts),uh1(Nvorts+1:2*Nvorts),gvals,gam,rval,pval);
k2 = dt*force_terms_multipole(mu,gam,rval,uh1,gvals,Nvorts,pval,tree_vals);
uh2 = u + k2/2;

tree_vals = multi_pole_kernel_build(uh2(1:Nvorts),uh2(Nvorts+1:2*Nvorts),gvals,gam,rval,pval);
k3 = dt*force_terms_multipole(mu,gam,rval,uh2,gvals,Nvorts,pval,tree_vals);
uh3 = u+k3;

tree_vals = multi_pole_kernel_build(uh3(1:Nvorts),uh3(Nvorts+1:2*Nvorts),gvals,gam,rval,pval);
k4 = dt*force_terms_multipole(mu,gam,rval,uh3,gvals,Nvorts,pval,tree_vals);

u = u + (k1 + 2*(k2+k3) + k4)/6;