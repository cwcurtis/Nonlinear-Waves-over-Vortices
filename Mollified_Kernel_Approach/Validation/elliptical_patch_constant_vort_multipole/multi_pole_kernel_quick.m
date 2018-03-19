function Kvec = multi_pole_kernel_quick(xpos,zpos,gvals,ep,pval,tree_val)

% Build kd-tree structure from xpos,zpos

Nvorts = length(xpos);

Kvec = zeros(Nvorts,2);
rcnt = 4;
ccnt = 4;
nblcks = rcnt*ccnt;

dx = tree_val{1,1}.dx;
dz = tree_val{1,1}.dz;

ctf = dx^2 + dz^2;

for jj=1:nblcks
    lnode = tree_val{jj,1};
    %linds = lnode.loc_list;
    linds = lnode.num_list;
    xloc = lnode.xpos;
    zloc = lnode.zpos;
    gloc = lnode.gvals;
    
    xcc = lnode.center;
    npts = lnode.tpts;
    cinds = [1:jj-1 jj+1:nblcks];
    Kfar = zeros(npts,2);
    %finds = zeros(Nvorts,1);
    finds = [];
    cmplst = [];
    for ll = 1:nblcks-1        
        cnode = tree_val{cinds(ll),1};        
        xcf = cnode.center;        
        nm_chldrn = cnode.no_chldrn;
        if cnode.tpts>0             
            dst = (xcc(1)-xcf(1))^2 + (xcc(2)-xcf(2))^2;
            if dst > ctf        
                kcur = cnode.kvals;
                rloc = -1./((xloc-xcf(1))+1i*(zloc-xcf(2)));
                qf = kcur(pval+1);
                for mm=1:pval
                    qf = kcur(pval+1-mm) + qf.*rloc;                        
                end
                qf = -rloc.*qf;
                Kfar = Kfar + [imag(qf) real(qf)];       
            elseif nm_chldrn > 0    
                % build nearest-neighbor list
                cmplst = [cmplst cinds(ll)];
            else
                % if we do not compute over something, then we save it for
                % descent.  
                %finds = finds + cnode.loc_list;             
                finds = [finds;cnode.num_list];
            end
        end
    end
    
    cmpnum = 1+length(cmplst);   
    if lnode.no_chldrn > 0
        dscnt_tree = cell(cmpnum,4);
        dscnt_tree(1,:) = tree_val(jj,2:5);
        dscnt_tree(2:cmpnum,:) = tree_val(cmplst,2:5);
        tvec = tree_traverser_quick(xpos,zpos,gvals,ep,pval,dscnt_tree,cmpnum,Nvorts,linds,finds);                 
    else
       %nninds = sum(tree_val(cmplst,1).loc_list);
       %nninds = logical(nninds+finds);
       nninds = [vertcat(tree_val(cmplst,1).num_list);finds];
       
       xlist = xpos(nninds);
       zlist = zpos(nninds);
       glist = gvals(nninds);
       
       tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,ep);         
    end    
    Kvec(linds,:) = Kfar + tvec;   
end