function Kvec = multi_pole_kernel_quick(xpos,zpos,gvals,gam,ep,pval,tree_val)

% Build kd-tree structure from xpos,zpos

Nvorts = length(xpos);
mlvl = floor(log(Nvorts)/log(4));

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

ctf = 2*sqrt(dx^2+dz^2);

clvl = 1;

for jj=1:nblcks
    lnode = tree_val{jj,1};
    linds = lnode.loc_list;
    iinds = lnode.int_list;        
    
    xloc = xpos(linds);
    zloc = zpos(linds);
    xcc = lnode.center;
    npts = lnode.tpts;
    cinds = [1:jj-1 jj+1:nblcks];
    Kfar = zeros(npts,2);
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
        end        
    end
    
    if(npts>mlvl)
       tvec = tree_traverser_quick(xpos,zpos,mlvl,clvl+1,gvals,gam,ep,pval,tree_val); 
    else
       nninds = logical(iinds-linds);
       xlist = xpos(nninds);
       zlist = zpos(nninds);
       gloc = gvals(linds);
       glist = gpos(nninds);
       tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,gam,ep); 
    end
    Kvec(linds,:) = Kfar+tvec;
end