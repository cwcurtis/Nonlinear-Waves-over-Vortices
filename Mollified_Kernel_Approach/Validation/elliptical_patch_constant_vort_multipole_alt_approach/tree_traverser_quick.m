function Kret = tree_traverser_quick(xpos,zpos,gvals,ep,pval,inc_tree,Nvorts,pinds)

% xpos - Global Nvorts x 1 array of x-positions.

% zpos - Global Nvorts x 1 array of z-positions.

% gvals - Global Nvorts x 1 array of circulation values.

% ep - mollification kernel radius

% pval - order of multipole expansion

% Nvorts - global number of vortices.  Provided in each function call so as
% to avoid unnecessary 'length' calls.  

% pinds - parent indices.

Kvec = zeros(Nvorts,2);

% Iterate over children, compute FMM approximation, and then descend or 
% finally compute where necessary.

for ll=1:4
    lleaf = inc_tree{ll};
    lnode = lleaf{1};
    npts = lnode.tpts;
        
    if npts>0 
        linds = lnode.num_list;
        xloc = lnode.xpos;
        zloc = lnode.zpos;
        gloc = lnode.gvals;
            
        if ~isempty(lnode.farlst)
            % This is a *highly* vectorized implementation of the FMM
            % approxmation.  
            kcurs = lnode.kcursf;
            xcfs = lnode.xcfs;
            no_toofar = size(xcfs,1);
            xlcmp = repmat(xloc+1i*zloc,1,no_toofar);
            xcfscmp = repmat((xcfs(:,1)+1i*xcfs(:,2)).',npts,1);
            rloc = -1./(xlcmp - xcfscmp);        
            rlocc = -1./(xlcmp - conj(xcfscmp));        
            kmat = repmat(kcurs(pval+1,:),npts,1); 
            qf = kmat;
            qfc = conj(kmat);
            for mm=1:pval
                kmat = repmat(kcurs(pval+1-mm,:),npts,1); 
                qf = kmat + qf.*rloc;                        
                qfc = conj(kmat) + qfc.*rlocc;                        
            end
            qf = sum(-rloc.*qf+rlocc.*qfc,2);
        else
            qf = zeros(npts,1);
        end
        
        if lnode.no_chldrn > 0
            tvec = tree_traverser_quick(xpos,zpos,gvals,ep,pval,lleaf(2:5),Nvorts,linds);            
        else
            % Here we finally compute at a terminal node.  
            nninds = lnode.nodscndlst;
            xlist = xpos(nninds);
            zlist = zpos(nninds);
            glist = gvals(nninds);            
            tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,ep);             
        end      
        Kvec(linds,:) = [imag(qf) real(qf)] + tvec;              
    end     
end
Kret = Kvec(pinds,:);