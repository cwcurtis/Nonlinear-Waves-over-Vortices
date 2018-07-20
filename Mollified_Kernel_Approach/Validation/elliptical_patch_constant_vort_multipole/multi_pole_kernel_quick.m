function Kvec = multi_pole_kernel_quick(xpos,zpos,gvals,ep,pval,tree_val)

% Build kd-tree structure from xpos,zpos

Nvorts = length(xpos);

Kvec = zeros(Nvorts,2);
rcnt = 4;
ccnt = 4;
nblcks = rcnt*ccnt;

for jj=1:nblcks
    lnode = tree_val{jj,1};
    npts = lnode.tpts;
    if npts > 0
        
        linds = lnode.num_list;
        xloc = lnode.xpos;
        zloc = lnode.zpos;
        gloc = lnode.gvals;
        
        if ~isempty(lnode.farlst)
            xcfs = lnode.xcfs;
            kcurs = lnode.kcursf;
            no_toofar = size(xcfs,1);
            xlcmp = repmat(xloc+1i*zloc,1,no_toofar);
            xcfscmp = repmat((xcfs(:,1)+1i*xcfs(:,2)).',npts,1);
            rloc = -1./(xlcmp - xcfscmp);        
            %rlocc = -1./(xlcmp - conj(xcfscmp));        
            kmat = repmat(kcurs(pval+1,:),npts,1); 
            qf = kmat;
            %qfc = conj(kmat);
            for mm=1:pval
                kmat = repmat(kcurs(pval+1-mm,:),npts,1); 
                qf = kmat + qf.*rloc;                        
                %qfc = conj(kmat) + qfc.*rlocc;                        
            end
            %qf = sum(-rloc.*qf+rlocc.*qfc,2);
            qf = sum(-rloc.*qf,2);
        else
            qf = zeros(npts,1);
        end
                
        if lnode.no_chldrn > 0
            tvec = tree_traverser_quick(xpos,zpos,gvals,ep,pval,tree_val(jj,2:5),Nvorts,linds);                 
        else
            nninds = lnode.nodscndlst;
            xlist = xpos(nninds);
            zlist = zpos(nninds);
            glist = gvals(nninds);
       
            tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,ep);              
        end
        Kvec(linds,:) = [imag(qf) real(qf)] + tvec;   
    end    
end