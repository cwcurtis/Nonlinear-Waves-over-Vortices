function Kchild = tree_traverser_list_maker(inc_prnt,inc_tree,pval,Nvorts,prvinds)

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

% Build information and interactions for and between each child.
flat_prnt = [inc_prnt{:}].';
lnodes = flat_prnt(:,1);
lheaders = [lnodes{:}];
lkids = flat_prnt(:,2:5);
has_kids = [lheaders.no_chldrn]>0;
nnlsts = zeros(4);

if ~isempty(prvinds)
    nodscndlst = {prvinds,prvinds,prvinds,prvinds};
else
    nodscndlst = cell(4,1);
end
            
for ll=1:4
    lnode = lheaders(ll);
    vcmp = [1:ll-1 ll+1:4];
    nnlsts(ll,vcmp) = has_kids(vcmp);
    numinds = lnode.num_list;
    if has_kids(ll) == 0 && lnode.tpts>0
       for jj=1:3
           nodscndlst{vcmp(jj)} = vertcat(nodscndlst{vcmp(jj)},numinds);       
       end
   end    
end

Kstore = cell(5,4);    

flat_tree = [inc_tree{:,:}].';
flag = 0;
if ~isempty(flat_tree)
    flag = 1;
    headers = flat_tree(:,1);
    headers = [headers{:}];
    cinds = [headers.tpts] > 0;
    indtmp = sum(cinds);
    inds_inc = (1:indtmp)';
    centers = [headers(cinds).center]';
    kvsary = [headers(cinds).kvals];
    cdcells = flat_tree(cinds,2:5);
    chldary = [headers(cinds).no_chldrn]' > 0;
    indtst = {headers.num_list};    
end

if flag == 1
    dx = lheaders(1).dx;
    dz = lheaders(1).dz; 
    ctf = dx^2+dz^2;
    nncntrs = [lheaders.center]';
    dxcntrs = repmat(nncntrs(:,1)',indtmp,1)-repmat(centers(:,1),1,4);
    dzcntrs = repmat(nncntrs(:,2)',indtmp,1)-repmat(centers(:,2),1,4);
    dsts = dxcntrs.^2 + dzcntrs.^2;
    toofar = dsts > ctf;
    tooclose = logical(1-toofar);
end

for ll=1:4
    lnode = lheaders(ll);    
    if lheaders(ll).tpts>0 && flag == 1        
        if ~isempty(toofar) 
            myfar = toofar(:,ll);
            lnode.farlst = inds_inc(myfar);
            lnode.kcursf = kvsary(:,myfar);
            lnode.xcfs = centers(myfar,:);
        end        
        myclose = tooclose(:,ll);
        cmplst = logical(myclose.*chldary);
        lnode.nearlst = inds_inc(cmplst);        
        frlst = logical(myclose - cmplst);        
        if sum(frlst)>0
           nodscndlst{ll} = vertcat(nodscndlst{ll},vertcat(indtst{frlst}));                
        end                        
    end
    
    if lnode.no_chldrn > 0
       nnnum = sum(nnlsts(ll,:));
       if flag == 1 && nnnum > 0
          cfnum = sum(cmplst);         
          cmpnum = cfnum+nnnum;            
          dscnt_tree = cell(cmpnum,4);            
          dscnt_tree(1:cfnum,:) = cdcells(cmplst,:);          
          dscnt_tree(cfnum+1:cmpnum,:) = lkids(logical(nnlsts(ll,:)),:);          
       elseif nnnum > 0
          dscnt_tree = lkids(logical(nnlsts(ll,:)),:);            
       else
          dscnt_tree = cell(1,4);
       end
       prnode = lkids(ll,:);      
       Kstore(2:5,ll) = tree_traverser_list_maker(prnode,dscnt_tree,pval,Nvorts,nodscndlst{ll});                                    
    else
       lnode.nodscndlst = nodscndlst{ll};            
    end
    Kstore{1,ll} = lnode;             
end

Kchild = {Kstore(:,1),Kstore(:,2),Kstore(:,3),Kstore(:,4)};
