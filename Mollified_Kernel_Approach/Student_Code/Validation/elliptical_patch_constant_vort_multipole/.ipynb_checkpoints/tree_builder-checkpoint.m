function Kloc = tree_builder(xpos,zpos,mlvl,clvl,nind,ctgry,gvals,xbnds,zbnds,gam,ep,pval,pnum,ll)

xmin = xbnds(1);
xmax = xbnds(2);
zmin = zbnds(1);
zmax = zbnds(2);
[ccnt,rcnt,nblcks,jvals] = top_props(ctgry);
dx = (xmax-xmin)/ccnt;
dz = (zmax-zmin)/rcnt;
       
jj = jvals(ll);
        
loc_data.bnum = nind + ll;
             
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
loc_data.bnum = nind + ll;
loc_data.parent = pnum;

if npts>0
   kvals = far_panel_comp(xloc,zloc,gloc,gam,xc,zc,pval);
else
   kvals = [];
end

loc_data.center = [xc;zc];
loc_data.kvals = kvals;
   
if npts > mlvl         
   nindn = nind + 16^clvl + 4*(ll-1);   
   Kloc = cell(5,1);
   Kloc{1} = loc_data;
   for ll=1:4
        Kchild = tree_builder(xpos,zpos,mlvl,clvl+1,nindn,ctgry,gvals,xbnds,zbnds,gam,ep,pval,nind+ll,ll);  
        Kloc{ll+1} = Kchild;
   end
else
   Kloc = {loc_data};
end
        
