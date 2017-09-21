function Kvec = multi_pole_kernel_one_cut(xpos,zpos,gvals,ep,Ntrunc)

% Build kd-tree structure from xpos,zpos
Nvorts = length(xpos);
mlvl = floor(log2(Nvorts));
Kx = zeros(Nvorts,1);
Kz = zeros(Nvorts,1);
Kvec = [Kx;Kz];
lvl = 2;

% Cut once to begin the process.  

xmin = min(xpos);
xmax = max(xpos);
zmin = min(zpos);
zmax = max(zpos);

xhf = (xmin+xmax)/2;
zhf = (zmin+zmax)/2;

ls = xpos < xhf;
rs = 1 - ls;
bs = zpos < zhf;
ts = 1 - bs;

