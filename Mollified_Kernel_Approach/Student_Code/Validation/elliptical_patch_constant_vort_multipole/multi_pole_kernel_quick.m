function Kvec = multi_pole_kernel_quick(xpos,zpos,gvals,gam,ep,pval,tree_val)

% Build kd-tree structure from xpos,zpos

Nvorts = length(xpos);

Kvec = zeros(Nvorts,2);

xmin = min(xpos);
xmax = max(xpos);
zmin = min(zpos);
zmax = max(zpos);

rcnt = 4;
ccnt = 4;
nblcks = rcnt*ccnt;

dx = (xmax-xmin)/ccnt;
dz = (zmax-zmin)/rcnt;

ctf = sqrt(dx^2+dz^2);

for jj=1:nblcks
    lnode = tree_val{jj,1};
    linds = lnode.loc_list;
    
    xloc = xpos(linds);
    zloc = zpos(linds);
    xcc = lnode.center;
    npts = lnode.tpts;
    cinds = [1:jj-1 jj+1:nblcks];
    Kfar = zeros(npts,2);
    cmplst = [];
    for ll = 1:nblcks-1        
        cnode = tree_val{cinds(ll),1};        
        xcf = cnode.center;        
        if norm(xcc-xcf) > ctf
            kcur = cnode.kvals;
            if ~isempty(kcur)
                q0 = kcur(1);                                
                zcn = (xloc-xcf(1))+1i*gam*(zloc-xcf(2));
                rloc = 1./zcn;
                zmat = zeros(npts,pval+1);
                zmat(:,1) = rloc;
                azcnsq = abs(rloc).^2;                
                for kk=2:pval+1
                    zmat(:,kk) = zmat(:,kk-1).*rloc; 
                end                
                qf = zmat*kcur;                
                Kfarn1 = q0*imag(zcn).*azcnsq + imag(qf);
                Kfarn2 = -q0*real(zcn).*azcnsq + real(qf);                
                Kfar = Kfar + [Kfarn1 Kfarn2];             
            end            
        else
            % build nearest-neighbor list
            cmplst = [cmplst cinds(ll)];
        end
    end
    
    cmpnum = length(cmplst);   
        
    if lnode.no_chldrn > 0
        cmpnum = length(cmplst);   
        dscnt_tree = cell(cmpnum+1,1);
        dscnt_tree{1} = {tree_val{1,2:5}};
        for mm=1:cmpnum
           dscnt_tree{mm+1} = {tree_val{cmplst(mm),2:5}}; 
        end
        tvec = tree_traverser_quick(xpos,zpos,gvals,dx/2,dz/2,gam,ep,pval,dscnt_tree,npts,cmpnum); 
    else
       nninds = zeros(Nvorts,1);
       for mm=1:cmpnum
           nninds = nninds + tree_val{cmplst(mm),1}.loc_list;
       end
       nninds = logical(nninds);
       xlist = xpos(nninds);
       zlist = zpos(nninds);
       gloc = gvals(linds);
       glist = gpos(nninds);
       tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,gam,ep);  
    end
    Kvec(linds,:) = Kfar+tvec;
end