function Kmaster = multi_pole_kernel_build(xpos,zpos,gvals,gam,ep,pval)

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
clvl = 1;

% 1st entry is struct with all necessary in cell information.  
% Next four entries are for children. 
Kmaster = cell(nblcks,5); 

dx = (xmax-xmin)/ccnt;
dz = (zmax-zmin)/rcnt;

for jj=1:nblcks
    
    col = mod(jj-1,ccnt);
    row = (jj-1-col)/ccnt;
    [dshift,ushift,lshift,rshift,ctgry] = shift_finder(col+1,row+1,rcnt,ccnt);
    
    ccl = col-lshift;
    ccr = col+rshift;
    rrt = row-dshift;
    rrb = row+ushift;
        
    xl = xmin + ccl*dx;
    xr = xmin + (ccr+1)*dx;
    zt = zmax - rrt*dz;
    zb = zmax - (rrb+1)*dz;
    xbnds = [xl;xr];
    zbnds = [zb;zt];   
    xccl = xmin + col*dx;
    xccr = xmin + (col+1)*dx;
    zcct = zmax - row*dz;
    zccb = zmax - (row+1)*dz;
    
    indsl = logical((xpos>=xccl).*(xpos<=xccr).*(zpos<=zcct).*(zpos>=zccb));    
    indsi = logical((xpos>=xl).*(xpos<=xr).*(zpos<=zt).*(zpos>=zb));        
    
    xc = (xccl+xccr)/2;
    zc = (zccb+zcct)/2;

    xloc = xpos(indsl);
    zloc = zpos(indsl);
    gloc = gvals(indsl);

    npts = sum(indsl);

    loc_data.loc_list = indsl;
    loc_data.int_list = indsi;
    loc_data.bnum = jj;
    loc_data.tpts = npts;
    
    if npts>0
        kvals = far_panel_comp(xloc,zloc,gloc,gam,xc,zc,pval);
    else
        kvals = [];
    end

    loc_data.center = [xc;zc];
    loc_data.kvals = kvals;
   
    Kmaster{jj,1} = loc_data;

    if npts > mlvl   
       nind = 16 + 4*(jj-1);
       for ll=1:4
           Kchild = tree_builder(xpos,zpos,mlvl,clvl+1,nind,ctgry,gvals,xbnds,zbnds,gam,ep,pval,jj,ll);           
           Kmaster{jj,ll+1} = Kchild;
       end
    end
      
end   