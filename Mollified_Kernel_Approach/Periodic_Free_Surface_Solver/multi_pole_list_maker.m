function tree_val = multi_pole_list_maker(tree_val,pval,Nvorts)

% Build kd-tree structure from xpos,zpos

rcnt = 4;
ccnt = 4;
nblcks = rcnt*ccnt;

dx = tree_val{1,1}.dx;
dz = tree_val{1,1}.dz;

ctf = dx^2 + dz^2;

for jj=1:nblcks
    lnode = tree_val{jj,1};
    farlst = [];
    nearlst = [];
    nodscndlst = [];
    centers = [];
    kvsary = [];
    
    lnode.farlst = farlst;
    lnode.nearlst = nearlst;
    lnode.nodscnlst = nodscndlst;
    lnode.xcfs = centers;
    lnode.kcursf = kvsary;
    
    if lnode.tpts > 0
        xcc = lnode.center;
        cinds = [1:jj-1 jj+1:nblcks];
        
        for ll = 1:nblcks-1        
            cnode = tree_val{cinds(ll),1};        
            xcf = cnode.center;        
            nm_chldrn = cnode.no_chldrn;
            if cnode.tpts>0             
                dst = (xcc(1)-xcf(1))^2 + (xcc(2)-xcf(2))^2;
                if dst > ctf        
                    farlst = [farlst;cinds(ll)];
                    centers = [centers; xcf(1) xcf(2)];
                    kvsary = [kvsary cnode.kvals];
                elseif nm_chldrn > 0    
                    % build nearest-neighbor list
                    nearlst = [nearlst;cinds(ll)];
                else
                    % if we do not compute over something, then we save it for
                    % descent.  
                    nodscndlst = [nodscndlst;cnode.num_list];
                end
            end
        end
        
        lnode.farlst = farlst;
        lnode.nearlst = nearlst;
        lnode.kcursf = kvsary;
        lnode.xcfs = centers;
        
        cmpnum = 1+length(nearlst);   
        if lnode.no_chldrn > 0
            dscnt_tree = cell(cmpnum,4);
            dscnt_tree(1,:) = tree_val(jj,2:5);
            dscnt_tree(2:cmpnum,:) = tree_val(nearlst,2:5);         
            
            tree_val(jj,2:5) = tree_traverser_list_maker(dscnt_tree,cmpnum,pval,Nvorts,nodscndlst);                  
        else
            lnode.nodscndlst = nodscndlst;                    
        end
    end    
    tree_val{jj,1} = lnode;
            
end