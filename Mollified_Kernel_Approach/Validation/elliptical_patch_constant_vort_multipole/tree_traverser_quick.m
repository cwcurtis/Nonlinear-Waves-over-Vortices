function Kret = tree_traverser_quick(xpos,zpos,gvals,ep,pval,inc_tree,nblcks,Nvorts,pinds,prvinds)

% xpos - Global Nvorts x 1 array of x-positions.

% zpos - Global Nvorts x 1 array of z-positions.

% gvals - Global Nvorts x 1 array of circulation values.

% ep - mollification kernel radius

% pval - order of multipole expansion

% inc_tree - inc(coming)_tree containing header inc_tree{1,1:4}{1:5}
% representing sibling children from common parent and
% inc_tree{2:nblcks,1:4}{1:5} representing all uncalculated interactions.  

% Nvorts - global number of vortices.  Provided in each function call so as
% to avoid unnecessary 'length' calls.  

% pinds - parent indices.

% prvinds - all those indices not computed over at the previous level due
% to terminal nodes that were too close at a given level.  These nodes are
% computed only at the end as part of the final, exact computation.  

Kvec = zeros(Nvorts,2);

% Create seperate lists for headers and children that everyone else in the
% local child list sees. 

fnum = nblcks*4;

centers = zeros(fnum,2);
indclls = zeros(Nvorts,fnum);
chldary = ones(fnum,1);
kvsary = zeros(pval+1,fnum);
cdcells = cell(fnum,4);
indtmp = 0;

for kk=2:nblcks
    branch = inc_tree(kk,:);
    
    for jj=1:4
        cnode = branch{jj}{1};
        if  cnode.tpts > 0
            indtmp = indtmp+1;   
            centers(indtmp,:) = cnode.center;
            kvsary(:,indtmp) = cnode.kvals;
            indclls(:,indtmp) = cnode.loc_list;   
            if cnode.no_chldrn > 0 % Note, we implicitly are pruning empty cells as we go here.
               cdcells(indtmp,:) = branch{jj}(2:5);                        
            else
               chldary(indtmp) = 0;            
            end
        end                
    end            
end

% Build information and interactions for and between each child.

mnodes = inc_tree(1,:);
dx = mnodes{1}{1}.dx;
dz = mnodes{1}{1}.dz; 
ctf = dx^2+dz^2;

finds = repmat(prvinds,1,4);
nnlsts = repmat(1:4,4,1);
nnlsts = nnlsts - diag(diag(nnlsts));
nncntrs = zeros(4,2);
llsts = zeros(Nvorts,4);

for ll=1:4
    lleaf = mnodes{ll};
    lnode = lleaf{1};
    llsts(:,ll) = lnode.loc_list;
    nncntrs(ll,:) = lnode.center;
    vcmp = [1:ll-1 ll+1:4];
    if lnode.no_chldrn == 0
       nnlsts(:,ll) = 0;
       finds(:,vcmp) = finds(:,vcmp) + repmat(llsts(:,ll),1,3);        
    end    
end

% Use very vectorized approach to determine distance of all child cells
% wiwth all potential interaction cells.  

dxcntrs = repmat(nncntrs(:,1)',indtmp,1)-repmat(centers(1:indtmp,1),1,4);
dzcntrs = repmat(nncntrs(:,2)',indtmp,1)-repmat(centers(1:indtmp,2),1,4);
dsts = dxcntrs.^2 + dzcntrs.^2;

toofar = dsts > ctf;
tooclose = logical(1-toofar);

% Iterate over children, compute FMM approximation, and then descend or 
% finally compute where necessary.

for ll=1:4
    lleaf = mnodes{ll};
    lnode = lleaf{1};
    npts = lnode.tpts;
        
    if npts>0 && ~isempty(toofar) 
        linds = logical(llsts(:,ll));
        
        % This is a *highly* vectorized implementation of the FMM
        % approxmation.  
        xloc = lnode.xpos;
        zloc = lnode.zpos;
        gloc = lnode.gvals;
        myfar = toofar(:,ll);
        myclose = tooclose(:,ll);
        kcurs = kvsary(:,myfar);
        xcfs = centers(myfar,:);
        no_toofar = size(xcfs,1);
        rloc = -1./(repmat(xloc+1i*zloc,1,no_toofar) - repmat((xcfs(:,1)+1i*xcfs(:,2)).',npts,1));        
        qf = repmat(kcurs(pval+1,:),npts,1);
        for mm=1:pval
            qf = repmat(kcurs(pval+1-mm,:),npts,1) + qf.*rloc;                        
        end
        qf = -sum(rloc.*qf,2);
        
        % look over nearest neighbors.
        nnlst = nnlsts(ll,:)~=0;                
        
        if lnode.no_chldrn > 0
            % Here we descend further down the tree.
            cmplst = logical(myclose.*chldary(1:indtmp));
            frlst = logical(myclose - cmplst);
        
            cfnum = sum(cmplst);
            nnnum = sum(nnlst);
            cmpnum = 1+cfnum+nnnum;
            dscnt_tree = cell(cmpnum,4);
            
            dscnt_tree(1,:) = mnodes{ll}(2:5);
            dscnt_tree(2:cfnum+1,:) = cdcells(cmplst,:);
            
            % Here we tack on nearest-neighbor among children interactions.
           
            if nnnum > 0
                nninds = 1:4;
                nninds = nninds(nnlst);            
                for mm=1:nnnum
                    dscnt_tree(cfnum+1+mm,:) = mnodes{nninds(mm)}(2:5);
                end
            end
            finds(:,ll) = finds(:,ll) + sum(indclls(:,frlst),2);
            tvec = tree_traverser_quick(xpos,zpos,gvals,ep,pval,dscnt_tree,cmpnum,Nvorts,linds,finds(:,ll));            
        else
            % Here we finally compute at a terminal node.  
            nninds = logical(sum(indclls(:,myclose),2)+sum(llsts(:,nnlst),2)+finds(:,ll));            
            xlist = xpos(nninds);
            zlist = zpos(nninds);
            glist = gvals(nninds);            
            tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,ep);             
        end      
        Kvec(linds,:) = [imag(qf) real(qf)] + tvec;                        
    end     
end
Kret = Kvec(pinds,:);