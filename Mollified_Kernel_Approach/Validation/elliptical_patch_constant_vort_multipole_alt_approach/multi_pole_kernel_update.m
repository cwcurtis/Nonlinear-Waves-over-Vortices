function tree_vals = multi_pole_kernel_update(u,gvals,pval,gam,Nvorts,tree_vals)

xpos = u(1:Nvorts);
zpos = gam*u(Nvorts+1:2*Nvorts);

ccnt = 4;
rcnt = 4;
nblcks = ccnt*rcnt;

mlvl = floor(log(Nvorts)/log(4));

% 1st entry is struct with all necessary in cell information.  
% Next four entries are for children. 

for jj=1:nblcks
    lnode = tree_vals{jj,1};
    numinds = lnode.num_list;
    
    xloc = xpos(numinds);
    zloc = zpos(numinds);
    gloc = gvals(numinds);
    
    lnode.xpos = xloc;
    lnode.zpos = zloc;
    lnode.gvals = gloc;
    npts = lnode.tpts;
    
    if npts>0
        lcenter = lnode.center;
        kvals = far_panel_comp(xloc,zloc,gloc,lcenter(1),lcenter(2),pval);
    else
        kvals = [];
    end

    lnode.kvals = kvals;
    tree_vals{jj,1} = lnode;
          
    if npts > mlvl   
       for ll=1:4
           dscnt_cell = tree_vals{jj,ll+1};
           Kchild = tree_traverser_update(xpos,zpos,mlvl,gvals,pval,dscnt_cell,Nvorts);           
           tree_vals{jj,ll+1} = Kchild;           
       end              
    end    
end   