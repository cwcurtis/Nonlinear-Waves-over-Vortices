function u = vort_update_multipole(mu,gam,rval,u,gvals,Nvorts,dt)

pval = 10;

% Build-out the tree
tree_vals = multi_pole_kernel_build(u,gvals,pval,gam,Nvorts);
tree_vals = multi_pole_list_maker(tree_vals,pval,Nvorts);
k1 = dt*force_terms_multipole(mu,gam,rval,u,gvals,Nvorts,pval,tree_vals);
uh1 = u + k1/2;

tree_vals = multi_pole_kernel_build(uh1,gvals,pval,gam,Nvorts);
tree_vals2 = multi_pole_list_maker(tree_vals,pval,Nvorts);
k2 = dt*force_terms_multipole(mu,gam,rval,uh1,gvals,Nvorts,pval,tree_vals2);
uh2 = u + k2/2;

tree_vals = multi_pole_kernel_build(uh2,gvals,pval,gam,Nvorts);
tree_vals3 = multi_pole_list_maker(tree_vals,pval,Nvorts);
k3 = dt*force_terms_multipole(mu,gam,rval,uh2,gvals,Nvorts,pval,tree_vals3);
uh3 = u + k3;

tree_vals = multi_pole_kernel_build(uh3,gvals,pval,gam,Nvorts);
tree_vals4 = multi_pole_list_maker(tree_vals,pval,Nvorts);
k4 = dt*force_terms_multipole(mu,gam,rval,uh3,gvals,Nvorts,pval,tree_vals4);

u = u + (k1 + 2*(k2+k3) + k4)/6;