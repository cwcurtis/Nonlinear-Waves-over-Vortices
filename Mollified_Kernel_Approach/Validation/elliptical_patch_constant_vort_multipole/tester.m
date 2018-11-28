function tester

npts = 10;
x = linspace(1,2,npts+1);
[xx,zz] = meshgrid(x,x);

nvorts = length(xx(:));
gvals = 1e-2*ones(nvorts,1);
pval = 1;
mx = 5;
ep = .01;

Kmaster = multi_pole_kernel_build(xx(:),zz(:),gvals,pval,mx,nvorts);
Kmaster = multi_pole_list_maker(Kmaster,pval,nvorts);
Kvec = multi_pole_kernel_quick(xx(:),zz(:),gvals,ep,pval,mx,Kmaster);
Kvecc = near_neighbor_comp_periodic(xx(:),zz(:),[],[],gvals,[],ep,mx);
disp(norm(Kvec-Kvecc))
disp([Kvec(:,1) Kvecc(:,1)])
disp([Kvec(:,2) Kvecc(:,2)])