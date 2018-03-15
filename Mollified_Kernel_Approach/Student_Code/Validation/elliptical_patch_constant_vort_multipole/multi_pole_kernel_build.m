function Kmaster = multi_pole_kernel_build(u,gvals,ep,pval,Nvorts)

xpos = u(1:Nvorts);
zpos = u(Nvorts+1:2*Nvorts);

% Build kd-tree structure from xpos,zpos
xmin = min(xpos);
xmax = max(xpos);
zmin = min(zpos);
zmax = max(zpos);

Nvorts = length(xpos);
mlvl = floor(log(Nvorts)/log(4));

xmin = xmin*(1 - sign(xmin)*.01);
xmax = xmax*(1 + sign(xmax)*.01);

zmin = zmin*(1 - sign(zmin)*.01);
zmax = zmax*(1 + sign(zmax)*.01);

ccnt = 4;
rcnt = 4;
nblcks = ccnt*rcnt;

% 1st entry is struct with all necessary in cell information.  
% Next four entries are for children. 
Kmaster = cell(nblcks,5); 

dx = (xmax-xmin)/ccnt;
dz = (zmax-zmin)/rcnt;

for jj=1:nblcks
    
    col = mod(jj-1,ccnt);
    row = (jj-1-col)/ccnt;
    
    xl = xmin + col*dx;
    xr = xmin + (col+1)*dx;
    zt = zmax - row*dz;
    zb = zmax - (row+1)*dz;
    
    indsl = logical((xpos>=xl).*(xpos<xr).*(zpos<=zt).*(zpos>zb));    
    
    xc = (xl+xr)/2;
    zc = (zb+zt)/2;

    xloc = xpos(indsl);
    zloc = zpos(indsl);
    gloc = gvals(indsl);

    npts = sum(indsl);
    
    loc_data.loc_list = indsl;
    loc_data.tpts = npts;
    loc_data.dx = dx;
    loc_data.dz = dz;
    
    if npts>0
        kvals = far_panel_comp(xloc,zloc,gloc,xc,zc,pval);
    else
        kvals = [];
    end

    loc_data.center = [xc;zc];
    loc_data.kvals = kvals;
   
    if npts > mlvl   
       loc_data.no_chldrn = 4;
       Kmaster{jj,1} = loc_data;
       dnx = dx/2;
       dnz = dz/2;
       for ll=1:4
           ncol = mod(ll-1,2);
           nrow = (ll-1-ncol)/2;
           xnl = xl + ncol*dnx;
           xnr = xl + (ncol+1)*dnx;
           znt = zt - nrow*dnz;
           znb = zt - (nrow+1)*dnz;
           Kchild = tree_builder(xpos,zpos,mlvl,gvals,xnl,xnr,znb,znt,ep,pval);           
           Kmaster{jj,ll+1} = Kchild;           
       end              
    else
       loc_data.no_chldrn = 0;
       Kmaster{jj,1} = loc_data; 
    end    
end   