function Kloc = tree_traverser_update(xpos,zpos,mlvl,gvals,pval,inc_cell,Nvorts)

lnode = inc_cell{1};

% 'numinds' is literal global position in the Nvorts long arrays
numinds = lnode.num_list;

% positions of points within this node.  
xloc = xpos(numinds);
zloc = zpos(numinds);
gloc = gvals(numinds);

npts = lnode.tpts;

% header data 
lnode.xpos = xloc;
lnode.zpos = zloc;
lnode.gvals = gloc;

if npts>0
   centers = lnode.center;
   kvals = far_panel_comp(xloc,zloc,gloc,centers(1),centers(2),pval);
else
   kvals = [];
end

lnode.kvals = kvals;

Kloc = cell(5,1);
Kloc{1} = lnode;
   
if npts > mlvl         
   for ll=1:4
       dscnt_cell = inc_cell{ll+1};
       Kchild = tree_traverser_update(xpos,zpos,mlvl,gvals,pval,dscnt_cell,Nvorts);  
       Kloc{ll+1} = Kchild;
   end
end
        
