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
    linds = lnode.loc_list;
    
    xloc = xpos(linds);
    zloc = zpos(linds);
    gloc = gvals(linds);
       
    xcc = lnode.center;
    npts = lnode.tpts;
    cinds = [1:jj-1 jj+1:nblcks];
    Kfar = zeros(npts,2);
    finds = zeros(Nvorts,1);
        
    cmplst = [];
    zmat = zeros(npts,pval+1);
    for ll = 1:nblcks-1        
        cnode = tree_val{cinds(ll),1};        
        xcf = cnode.center;        
        nm_chldrn = cnode.no_chldrn;
        if cnode.tpts>0             
            dst = (xcc(1)-xcf(1))^2 + (xcc(2)-xcf(2))^2;
            if dst > ctf        
                kcur = cnode.kvals;
                zcn = (xloc-xcf(1))+1i*(zloc-xcf(2));
                rloc = -1./zcn;
                               
                zmat(:,1) = rloc;
                for kk=2:pval+1
                    zmat(:,kk) = zmat(:,kk-1).*rloc; 
                end                
                qf = zmat*kcur;                
                
                Kfar = Kfar - [imag(qf) real(qf)];       
            elseif nm_chldrn > 0    
                % build nearest-neighbor list
                cmplst = [cmplst cinds(ll)];
            else
                finds = finds + cnode.loc_list;                 
            end
        end
    end
    
    % compute over anything we do not descend on 
    if sum(finds)>0
        finds = logical(finds);
        Kndscnd = far_panel_exact_comp(xloc,zloc,xpos(finds),zpos(finds),gvals(finds),ep);              
        Kfar = Kfar + Kndscnd;
    end
        
    cmpnum = length(cmplst);   
    cmpinds = zeros(Nvorts,1);
    if lnode.no_chldrn > 0
        dscnt_tree = cell(cmpnum+1,1);
        dscnt_tree{1} = {tree_val{jj,2:5}};
        for mm=1:cmpnum
           dscnt_tree{mm+1} = {tree_val{cmplst(mm),2:5}};            
           cmpinds = cmpinds + tree_val{cmplst(mm),1}.loc_list;           
        end
        tvec = tree_traverser_quick(xpos,zpos,gvals,ep,pval,dscnt_tree,cmpnum,Nvorts,linds);                 
    else
       nninds = zeros(Nvorts,1);
       for mm=1:cmpnum
           nninds = nninds + tree_val{cmplst(mm),1}.loc_list;
       end
       nninds = logical(nninds);
       
       xlist = xpos(nninds);
       zlist = zpos(nninds);
       glist = gvals(nninds);
       
       tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,ep);         
    end    
    Kvec(linds,:) = Kfar + tvec;   
end