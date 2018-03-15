function Kloc = tree_builder(xpos,zpos,mlvl,gvals,xl,xr,zb,zt,ep,pval)

dx = xr-xl;
dz = zt-zb;

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

Kloc = cell(5,1);

if npts > mlvl         
   loc_data.no_chldrn = 4;
   Kloc{1} = loc_data;
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
       Kloc{ll+1} = Kchild;
   end
else
   loc_data.no_chldrn=0;
   Kloc{1} = loc_data;
end
        
