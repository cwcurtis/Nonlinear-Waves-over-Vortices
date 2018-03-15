function Kret = tree_traverser_quick(xpos,zpos,gvals,ep,pval,inc_tree,nblcks,Nvorts,pinds)

Kvec = zeros(Nvorts,2);

%Create seperate lists for headers and children that everyone else in the
%local list sees. 
fnum = nblcks*4;

tptsary = zeros(fnum,1);
centers = zeros(fnum,2);
indclls = cell(fnum,1);
chldary = zeros(fnum,1);
cdcells = cell(fnum,1);
kvsary = cell(fnum,1);

for kk=2:nblcks+1
    for jj=1:4
        cnode = inc_tree{kk}{jj}{1};
        ctpts = cnode.tpts;
        tptsary((kk-2)*4 + jj) = cnode.tpts;
        
        if ctpts > 0
            chldary((kk-2)*4 + jj) = cnode.no_chldrn;
            indclls{(kk-2)*4 + jj} = cnode.loc_list;        
            centers((kk-2)*4 + jj,:) = cnode.center;
            kvsary{(kk-2)*4 + jj} = cnode.kvals;
        end
        if cnode.no_chldrn > 0
            cdcells{(kk-2)*4 + jj} = {inc_tree{kk}{jj}{2:5}};
        end
    end            
end
    
for ll=1:4
    lleaf = inc_tree{1}{ll};
    lnode = lleaf{1};
    npts = lnode.tpts;
    vals = (1:4)';
    vcmp = vals([1:ll-1 ll+1:4]);
    if npts>0
        dx = lnode.dx;
        dz = lnode.dz;
        ctf = dx^2+dz^2;
        linds = lnode.loc_list;
        xloc = xpos(linds);
        zloc = zpos(linds);
        gloc = gvals(linds);
        finds = zeros(Nvorts,1);
        xcc = lnode.center;
        Kfar = zeros(npts,2);
        cmplst = [];
        
        for kk=1:4*nblcks
            if tptsary(kk)>0
               xcf = centers(kk,:);
               nm_chldrn = chldary(kk);
               dst = (xcc(1)-xcf(1))^2 + (xcc(2)-xcf(2))^2;
               if  dst > ctf 
                   kcur = kvsary{kk};                   
                   zcn = (xloc-xcf(1))+1i*(zloc-xcf(2));
                   rloc = -1./zcn;
                   
                   qf = kcur(pval+1);
                   for mm=1:pval
                       qf = kcur(pval+1-mm) + qf.*rloc;                        
                   end
                   qf = rloc.*qf;
                   qfi = -imag(qf);
                   qfr = -real(qf);
                      
                   % Maybe someday this mex-file will be faster.
                   %{ 
                   rlr = real(rloc);
                   rli = imag(rloc);
                   kcr = -real(kcur);
                   kci = -imag(kcur);
                   
                   [qfr,qfi] = qf_comp(rlr,rli,kcr,kci);
                   %}
                   
                   Kfar(:,1) = Kfar(:,1) + qfi;
                   Kfar(:,2) = Kfar(:,2) + qfr;
               elseif nm_chldrn > 0
                   % build nearest-neighbor list
                   cmplst = [cmplst kk];
               else
                   finds = finds + indclls{kk};                   
               end
            end                      
        end
        
        % look over nearest neighbors.
        nnlst = [];
        for jj=1:3
           cnode = inc_tree{1}{vcmp(jj)}{1};
           if cnode.no_chldrn > 0 
               nnlst = [nnlst vcmp(jj)];               
           else
               finds = finds + cnode.loc_list;               
           end
        end
        
        % compute over all those things too close for multipole, and with
        % no children to descend over.  
        if sum(finds)>0
            finds = logical(finds);
            Kndscnd = far_panel_exact_comp(xloc,zloc,xpos(finds),zpos(finds),gvals(finds),ep);              
            Kfar = Kfar + Kndscnd;
        end
        
        cfnum = length(cmplst);
        cmpnum = cfnum+length(nnlst);
        if lnode.no_chldrn > 0
            dscnt_tree = cell(cmpnum+1,1);
            dscnt_tree{1} = {lleaf{2:5}};                
            for mm=1:cmpnum
                if mm<=cfnum
                    dscnt_tree{mm+1} = cdcells{cmplst(mm)};
                else
                    mms = mm-cfnum;
                    dscnt_tree{mm+1} = {inc_tree{1}{nnlst(mms)}{2:5}};
                end
            end                        
            tvec = tree_traverser_quick(xpos,zpos,gvals,ep,pval,dscnt_tree,cmpnum,Nvorts,linds);            
        else
            nninds = zeros(Nvorts,1);
            for mm=1:cmpnum
                if mm<=cfnum
                    nninds = nninds + indclls{cmplst(mm)};                
                else
                    mms = mm-cfnum;
                    nninds = nninds + inc_tree{1}{nnlst(mms)}{1}.loc_list;
                end                
            end
          
            nninds = logical(nninds);
            
            xlist = xpos(nninds);
            zlist = zpos(nninds);
            glist = gvals(nninds);
            
            tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,ep);             
        end      
        Kvec(linds,:) = Kfar + tvec;                        
    end     
end
Kret = Kvec(pinds,:);