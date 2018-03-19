function Kloc = tree_builder(xpos,zpos,mlvl,gvals,xl,xr,zb,zt,ep,pval,pinds,Nvorts)

dx = xr-xl;
dz = zt-zb;

%So that we do not have to keep searching over Nvorts long lists, we pass
%parent indices in.  
xcur = xpos(pinds);
zcur = zpos(pinds);

inds_in_cell = logical((xcur>=xl).*(xcur<xr).*(zcur<=zt).*(zcur>zb));

xc = (xl+xr)/2;
zc = (zb+zt)/2;

% 'numinds' is literal global position in the Nvorts long arrays
numinds = pinds(inds_in_cell);

% 'indsl' is a logical indexing Nvorts long variant of 'numinds' so that we
% can simply sum positions if we like.  
indsl = zeros(Nvorts,1);
indsl(numinds) = 1;

% positions of points within this node.  
xloc = xpos(numinds);
zloc = zpos(numinds);
gloc = gvals(numinds);

npts = sum(indsl);

% header data 
loc_data.loc_list = indsl;
loc_data.num_list = pinds(inds_in_cell);
loc_data.xpos = xloc;
loc_data.zpos = zloc;
loc_data.gvals = gloc;
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
       Kchild = tree_builder(xpos,zpos,mlvl,gvals,xnl,xnr,znb,znt,ep,pval,numinds,Nvorts);  
       Kloc{ll+1} = Kchild;
   end
else
   loc_data.no_chldrn=0;
   Kloc{1} = loc_data;
end
        
