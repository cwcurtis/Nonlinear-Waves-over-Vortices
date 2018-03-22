function Kchild = tree_traverser_list_maker(inc_tree,nblcks,pval,Nvorts,prvinds)

% inc_tree - inc(coming)_tree containing header inc_tree{1,1:4}{1:5}
% representing sibling children from common parent and
% inc_tree{2:nblcks,1:4}{1:5} representing all uncalculated interactions.  

% Nvorts - global number of vortices.  Provided in each function call so as
% to avoid unnecessary 'length' calls.  

% prvinds - all those indices not computed over at the previous level due
% to terminal nodes that were too close at a given level.  These nodes are
% computed only at the end as part of the final, exact computation.  

% Create seperate lists for headers and children that everyone else in the
% local child list sees. 

Kstore = cell(5,4);    
fnum = nblcks*4;

centers = zeros(fnum,2);
indtst = cell(fnum,1);
chldary = ones(fnum,1);
cdcells = cell(fnum,4);
kvsary = zeros(pval+1,fnum);

indtmp = 0;

for kk=2:nblcks
    branch = inc_tree(kk,:);
    for jj=1:4
        cnode = branch{jj}{1};        
        if  cnode.tpts > 0
            indtmp = indtmp+1;   
            centers(indtmp,:) = cnode.center;
            indtst{indtmp} = cnode.num_list;
            kvsary(:,indtmp) = cnode.kvals;            
            if cnode.no_chldrn > 0 
               cdcells(indtmp,:) = branch{jj}(2:5);                                       
            else
               chldary(indtmp) = 0;            
            end
        end        
    end            
end

inds_inc = (1:indtmp)';

% Build information and interactions for and between each child.

mnodes = inc_tree(1,:);
dx = mnodes{1}{1}.dx;
dz = mnodes{1}{1}.dz; 
ctf = dx^2+dz^2;

if ~isempty(prvinds)
    nodscndlst = {prvinds,prvinds,prvinds,prvinds};
else
    nodscndlst = cell(4,1);
end

nnlsts = repmat(1:4,4,1);
nnlsts = nnlsts - diag(diag(nnlsts));
nncntrs = zeros(4,2);
llsts = cell(4,1);
nninds = (1:4)';
            
for ll=1:4
    lnode = mnodes{ll}{1};
    numinds = lnode.num_list;
    llsts{ll} = numinds;
    nncntrs(ll,:) = lnode.center;
    vcmp = [1:ll-1 ll+1:4];
    if lnode.no_chldrn == 0 
       nnlsts(:,ll) = 0;
       if ~isempty(numinds)
            for jj=1:3
                nodscndlst{vcmp(jj)} = vertcat(nodscndlst{vcmp(jj)},numinds);
            end
       end
   end    
end

% Use very vectorized approach to determine distance of all child cells
% wiwth all potential interaction cells.  

dxcntrs = repmat(nncntrs(:,1)',indtmp,1)-repmat(centers(1:indtmp,1),1,4);
dzcntrs = repmat(nncntrs(:,2)',indtmp,1)-repmat(centers(1:indtmp,2),1,4);
dsts = dxcntrs.^2 + dzcntrs.^2;

toofar = dsts > ctf;
tooclose = logical(1-toofar);

for ll=1:4
    lnode = mnodes{ll}{1};    
    if lnode.tpts>0 
        
        if ~isempty(toofar) 
            myfar = toofar(:,ll);
            lnode.farlst = inds_inc(myfar);
            lnode.kcursf = kvsary(:,myfar);
            lnode.xcfs = centers(myfar,:);
        end
        
        myclose = tooclose(:,ll);
        cmplst = logical(myclose.*chldary(1:indtmp));
        lnode.nearlst = inds_inc(cmplst);
        
        frlst = logical(myclose - cmplst);        
        if sum(frlst)>0
           nodscndlst{ll} = vertcat(nodscndlst{ll},vertcat(indtst{frlst}));                
        end                
        
        if lnode.no_chldrn > 0
            % Here we descend further down the tree.
            % look over nearest neighbors.
            nnlst = nnlsts(ll,:)~=0;                
            cfnum = sum(cmplst);
            nnnum = sum(nnlst);
            cmpnum = 1+cfnum+nnnum;
            dscnt_tree = cell(cmpnum,4);            
            dscnt_tree(1,:) = mnodes{ll}(2:5);
            dscnt_tree(2:cfnum+1,:) = cdcells(cmplst,:);            
            % Here we tack on nearest-neighbor among children interactions.           
            if nnnum > 0
                nnindr = nninds(nnlst);            
                mnodesrem = mnodes(nnindr);            
                for mm=1:nnnum
                    dscnt_tree(cfnum+1+mm,:) = mnodesrem{mm}(2:5);
                end                
            end            
            Kstore(2:5,ll) = tree_traverser_list_maker(dscnt_tree,cmpnum,pval,Nvorts,nodscndlst{ll});                                    
        else
            lnode.nodscndlst = nodscndlst{ll};            
        end               
    end
    Kstore{1,ll} = lnode;             
end

Kchild = {Kstore(:,1),Kstore(:,2),Kstore(:,3),Kstore(:,4)};
