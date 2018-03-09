function Kvec = tree_traverser_quick(xpos,zpos,mlvl,clvl,gvals,gam,ep,pval,pnum,tree_val)



for ll=1:4
        
    jj = jvals(ll);
        
    Kloc{ll,1} = nind + ll;
             
    col = mod(jj-1,ccnt);
    row = (jj-1-col)/ccnt;
    [dshift,ushift,lshift,rshift,ctgry] = shift_finder(col+1,row+1,rcnt,ccnt);
    
    ccl = col-lshift;
    ccr = col+rshift;
    rrt = row-dshift;
    rrb = row+ushift;
           
    cmpvals = ones(rcnt,ccnt);
    cmpvals(rrt+1:rrb+1,ccl+1:ccr+1)=0;
    mask = cmpvals(:);
    inds = mask == 1;
    
    crvals = rvals(inds);
    ccvals = cvals(inds);
           
    xl = xmin + (col-lshift)*dx;
    xr = xmin + (col+rshift+1)*dx;
    zt = zmax - (row-dshift)*dz;
    zb = zmax - (row+ushift+1)*dz;
    xbndsloc = [xl;xr];
    zbndsloc = [zb;zt];            
    
    xccl = xmin + col*dx;
    xccr = xmin + (col+1)*dx;
    zcct = zmax - row*dz;
    zccb = zmax - (row+1)*dz;
    
    indsi = logical((xpos>=xl).*(xpos<=xr).*(zpos<=zt).*(zpos>=zb));
    indsl = logical((xpos>=xccl).*(xpos<=xccr).*(zpos<=zcct).*(zpos>=zccb));
    
    farinds = logical(1-indsi);    
        
    if sum(indsi)>0
    
        xcells = xpos(indsi);
        xfar = xpos(farinds);
        
        zcells = zpos(indsi);
        zfar = zpos(farinds);

        gcells = gvals(indsi);
        gfar = gvals(farinds);
        
        indsnn = logical(indsi-indsl);

        xloc = xpos(indsl);
        zloc = zpos(indsl);
        gloc = gvals(indsl);

        xlist = xpos(indsnn);
        zlist = zpos(indsnn);
        glist = gvals(indsnn);
    
        npartsloc = sum(indsl);       

        % just need to loop over the far-field here
        Kfar = zeros(npartsloc,2);
        tvec = zeros(npartsloc,2);
        if ~isempty(crvals)
            for mm=1:length(crvals)
                rnuml = crvals(mm);
                cnuml = ccvals(mm);      
                if isempty(Kcomp{rnuml,cnuml})                       
                   xll = xmin + dx*(cnuml-1);
                   xrl = xmin + dx*(cnuml);                
                   ztl = zmax - dz*(rnuml-1);
                   zbl = zmax - dz*(rnuml);                
                   indsf = logical((xfar>=xll).*(xfar<=xrl).*(zfar>=zbl).*(zfar<=ztl));        
                   xfarl = xfar(indsf);
                   zfarl = zfar(indsf);
                   gfarl = gfar(indsf);                        
                   indsc = logical(1-indsf);
                   xfar = xfar(indsc);
                   zfar = zfar(indsc);
                   gfar = gfar(indsc);
                   xc = (xll+xrl)/2;
                   zc = (ztl+zbl)/2;                           
                   if ~isempty(xfarl)
                       Kcomp{rnuml,cnuml} = far_panel_comp(xfarl,zfarl,gfarl,gam,xc,zc,pval);  
                       Cents{rnuml,cnuml} = [xc;zc];
                       Kloc{ll,2} = Kcomp{rnuml,cnuml};
                   end
                end
                
                if ~isempty(Kcomp{rnuml,cnuml})
                    vc = Cents{rnuml,cnuml};
                    qvalsl = Kcomp{rnuml,cnuml};
                    q0 = qvalsl(1);                
                
                    zcn = (xloc-vc(1))+1i*gam*(zloc-vc(2));
                    rloc = 1./zcn;
                    zmat = zeros(npartsloc,pval+1);
                    zmat(:,1) = rloc;
                    azcnsq = abs(rloc).^2;
                
                    for kk=2:pval+1
                        zmat(:,kk) = zmat(:,kk-1).*rloc; 
                    end
                
                    qf = zmat*qvalsl;
                
                    Kfarn1 = q0*imag(zcn).*azcnsq + imag(qf);
                    Kfarn2 = -q0*real(zcn).*azcnsq + real(qf);                
                    Kfar = Kfar + [Kfarn1 Kfarn2];             
                end
            end
        end
        
        if npartsloc > mlvl         
            nindn = nind + 16^clvl + 4*(ll-1);   
            [tvec,Kchild] = tree_traverser_quick(xcells,zcells,npartsloc,mlvl,clvl+1,nindn,ctgry,gcells,xbndsloc,zbndsloc,gam,ep,pval);  
            Kloc{ll,3} = Kchild;
        elseif npartsloc > 0            
            % Find nearest-neighbor interactions      
            tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,gam,ep);                    
        end
        Kvec(ind_cnt+1:ind_cnt+npartsloc,:) = Kfar + tvec;  
        ind_cnt = ind_cnt + npartsloc;    
    end       
end