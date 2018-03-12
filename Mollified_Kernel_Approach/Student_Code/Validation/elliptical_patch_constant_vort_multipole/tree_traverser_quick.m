function Kvec = tree_traverser_quick(xpos,zpos,gvals,dx,dz,gam,ep,pval,inc_tree,totpts,nblcks)

ctf = sqrt(dx^2+dz^2);
ind_cnt = 0;
Kvec = zeros(totpts,2);
Nvorts = length(xpos);

for ll=1:4
    lleaf = inc_tree{1}{ll};
    lnode = lleaf{1};
    linds = lnode.loc_list;
    
    xloc = xpos(linds);
    zloc = zpos(linds);
    
    xcc = lnode.center;
    npts = lnode.tpts;
    Kfar = zeros(npts,2);
    cmplst = [];
    
    for kk=2:nblcks
        for jj=1:4
            if ~isempty(inc_tree{kk}{jj})
                cnode = inc_tree{kk}{jj}{1};        
                if ~isempty(cnode.kvals)
                    xcf = cnode.center;        
                    if norm(xcc-xcf) > ctf
                       kcur = cnode.kvals;
                       q0 = kcur(1);                                
                       zcn = (xloc-xcf(1))+1i*gam*(zloc-xcf(2));
                       rloc = 1./zcn;
                       zmat = zeros(npts,pval+1);
                       zmat(:,1) = rloc;
                       azcnsq = abs(rloc).^2;                
                       for mm=2:pval+1
                           zmat(:,mm) = zmat(:,mm-1).*rloc; 
                       end                
                       qf = zmat*kcur;                
                       Kfarn1 = q0*imag(zcn).*azcnsq + imag(qf);
                       Kfarn2 = -q0*real(zcn).*azcnsq + real(qf);                
                       Kfar = Kfar + [Kfarn1 Kfarn2];             
                    else
                    % build nearest-neighbor list
                       cmplst = [cmplst [kk;jj]];
                    end
                end
            end
        end
    end
    
    cmpnum = size(cmplst,2);
    
    if lnode.no_chldrn > 0
       dscnt_tree = cell(cmpnum+1,1);
       dscnt_tree{1} = {lleaf{2:5}};
       for mm=1:cmpnum
           dscnt_tree{mm+1} = {inc_tree{cmplst(1,mm)}{cmplst(2,mm)}{2:5}};
       end
       tvec = tree_traverser_quick(xpos,zpos,gvals,dx/2,dz/2,gam,ep,pval,dscnt_tree,npts,cmpnum);
    else
       nninds = zeros(Nvorts,1);
       for mm=1:cmpnum
           nninds = nninds + inc_tree{cmplst(1,mm)}{cmplst(2,mm)}{1}.loc_list;
       end 
       nninds = logical(nninds);
       xlist = xpos(nninds);
       zlist = zpos(nninds);
       gloc = gvals(linds);
       glist = gvals(nninds);
       tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,gam,ep); 
    end
    
    Kvec(ind_cnt+1:ind_cnt+npts,:) = Kfar + tvec;  
    
    ind_cnt = ind_cnt + npts;            
end

disp([ind_cnt totpts])