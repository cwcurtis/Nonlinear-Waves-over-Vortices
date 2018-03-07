function Kvec = multi_pole_kernel(xpos,zpos,gvals,gam,ep,pval)

% Build kd-tree structure from xpos,zpos
Nvorts = length(xpos);
mlvl = floor(log(Nvorts)/log(4));

Kvec = zeros(Nvorts,2);
xmin = min(xpos);
xmax = max(xpos);
zmin = min(zpos);
zmax = max(zpos);

xmin = xmin*(1 - sign(xmin)*.01);
xmax = xmax*(1 + sign(xmax)*.01);

zmin = zmin*(1 - sign(zmin)*.01);
zmax = zmax*(1 + sign(zmax)*.01);

ccnt = 4;
rcnt = 4;
nblcks = 16;

rvals = kron(ones(ccnt,1),(1:rcnt)');
cvals = kron((1:ccnt)',ones(rcnt,1));

Kcomp = cell(mlvl,6,6);
Cents = cell(mlvl,6,6);
clvl = 1;

dx = (xmax-xmin)/ccnt;
dz = (zmax-zmin)/rcnt;

vvecs = zeros(Nvorts,ccnt);
hvecs = zeros(Nvorts,rcnt);
bvecs = zeros(Nvorts,nblcks);

ilists = zeros(Nvorts,nblcks);
farfield = zeros(Nvorts,nblcks);

for jj=1:ccnt
    xl = xmin + (jj-1)*dx;
    xr = xl + dx;
    if jj<ccnt
        vvecs(:,jj) = (xpos>=xl).*(xpos<xr);    
    else
        vvecs(:,jj) = (xpos>=xl).*(xpos<=xmax);    
    end
end

for jj=1:rcnt
    zt = zmax - (jj-1)*dz;
    zb = zt - dz;
    if jj<rcnt
        hvecs(:,jj) = (zpos>zb).*(zpos<=zt); 
    else
        hvecs(:,jj) = (zpos>=zmin).*(zpos<=zt); 
    end
end

for jj = 1:nblcks
    col = mod(jj-1,ccnt);
    row = (jj-1-col)/ccnt;
    bvecs(:,jj) = hvecs(:,row+1).*vvecs(:,col+1);    
end

for jj=1:nblcks
    
    col = mod(jj-1,ccnt);
    row = (jj-1-col)/ccnt;
    [dshift,ushift,lshift,rshift,ctgry] = shift_finder(col+1,row+1,rcnt,ccnt);
    ccl = col-lshift;
    ccr = col+rshift;
    rrt = row-dshift;
    rrb = row+ushift;
    
    %for rnum = rrt:rrb
    %    for cnum = ccl:ccr
    %        ilists(:,jj) = ilists(:,jj) + bvecs(:,ccnt*rnum+cnum+1);            
    %    end
    %end
    
    crinds = ccnt*(rrt:rrb)'*ones(1,length(ccl:ccr)) + ones(length(rrt:rrb),1)*(ccl:ccr) + ones(length(rrt:rrb),length(ccl:ccr));
    ilists(:,jj) = sum(bvecs(:,sort(crinds(:))),2);
        
    xl = xmin + (col-lshift)*dx;
    xr = xmin + (col+rshift+1)*dx;
    zt = zmax - (row-dshift)*dz;
    zb = zmax - (row+ushift+1)*dz;
    xbnds = [xl;xr];
    zbnds = [zb;zt];            
    
    % Find complement box coordinates so as to loop over far-field boxes.
    
    cmpvals = ones(rcnt,ccnt);
    cmpvals(rrt+1:rrb+1,ccl+1:ccr+1)=0;
    mask = cmpvals(:);
    inds = mask == 1;
    
    crvals = rvals(inds);
    ccvals = cvals(inds);
    
    nparts = sum(bvecs(:,jj));
    
    if nparts > 0
        farfield(:,jj) = 1 - ilists(:,jj);
    
        locinds = bvecs(:,jj)~=0;
        lisinds = ilists(:,jj)-bvecs(:,jj)~=0;
        tot_inds = ilists(:,jj)~=0;
        farinds = farfield(:,jj)~=0;
    
        xloc = xpos(locinds);
        xlist = xpos(lisinds);
        xcells = xpos(tot_inds);
        xfar = xpos(farinds);
    
        zloc = zpos(locinds);
        zlist = zpos(lisinds);
        zcells = zpos(tot_inds);
        zfar = zpos(farinds);
    
        gloc = gvals(locinds);
        glist = gvals(lisinds);
        gcells = gvals(tot_inds);
        gfar = gvals(farinds);    
    
        % just need to loop over the far-field here
        Kfar = zeros(length(xloc),2);
        if ~isempty(crvals)
            for mm=1:length(crvals)
                rnuml = crvals(mm);
                cnuml = ccvals(mm);      
                if isempty(Kcomp{clvl,rnuml,cnuml})   
                        xll = xmin + dx*(cnuml-1);
                        xrl = xmin + dx*(cnuml);                
                        ztl = zmax - dz*(rnuml-1);
                        zbl = zmax - dz*(rnuml);                
                        inds = logical((xfar>=xll).*(xfar<=xrl).*(zfar>=zbl).*(zfar<=ztl));                
                        xc = (xll+xrl)/2;
                        zc = (ztl+zbl)/2;                           
                        xfarl = xfar(inds);
                        zfarl = zfar(inds);
                        gfarl = gfar(inds);
                        Kcomp{clvl,rnuml,cnuml} = far_panel_comp(xfarl,zfarl,gfarl,gam,xc,zc,pval);  
                        Cents{clvl,rnuml,cnuml} = [xc;zc];
                end
                vc = Cents{clvl,rnuml,cnuml};
                zcn = xloc+1i*gam*zloc - (vc(1)+1i*gam*vc(2));
                rloc = 1./zcn;
                azcnsq = abs(zcn).^2;
                qvalsl = Kcomp{clvl,rnuml,cnuml};
                q0 = qvalsl(1);                
                qf = qvalsl(2)*rloc;
                for kk=3:length(qvalsl)
                    rloc = rloc./zcn;
                    qf = qf + qvalsl(kk)*rloc;
                end
                
                Kfarn1 = q0*imag(zcn)./azcnsq + imag(qf);
                Kfarn2 = -q0*real(zcn)./azcnsq + real(qf);                
                Kfar = Kfar + [Kfarn1 Kfarn2];  
            end
        end
        %Kfart = far_panel_exact_comp(xloc,zloc,xfar,zfar,gfar,gam);

        %plot(xloc,log10(abs(Kfar(:,1))),'k--',xloc,log10(abs(Kfart(:,1))),'b-','LineWidth',2)
        %pause
        
        if nparts > mlvl
            tvec = tree_traverser(xcells,zcells,nparts,mlvl,clvl+1,Kcomp,Cents,ctgry,gcells,xbnds,zbnds,gam,ep,pval);            
        else
            tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,gam,ep);            
        end
        Kvec(locinds,:) = Kfar + tvec;        
    end
end