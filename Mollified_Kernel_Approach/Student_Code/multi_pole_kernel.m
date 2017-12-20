function Kvec = multi_pole_kernel(xpos,zpos,gvals,Ntrunc,gam,ep)

% Build kd-tree structure from xpos,zpos
Nvorts = length(xpos);
mlvl = floor(log2(Nvorts));
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
    
    for rnum = row-dshift:row+ushift
        for cnum = col-lshift:col+rshift
            ilists(:,jj) = ilists(:,jj) + bvecs(:,ccnt*rnum+cnum+1);            
        end
    end
    
    xl = xmin + (col-lshift)*dx;
    xr = xmin + (col+rshift+1)*dx;
    zt = zmax - (row-dshift)*dz;
    zb = zmax - (row+ushift+1)*dz;
    xbnds = [xl;xr];
    zbnds = [zb;zt];            
    
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
    
        Kfar = far_panel_comp(xloc,zloc,xfar,zfar,gfar,gam,2);        
        %Kfar_ex = far_panel_exact_comp(xloc,zloc,xfar,zfar,gfar,Ntrunc,gam,ep);        
        %plot(1:nparts,Kfar(:,1),'--',1:nparts,Kfar_ex(:,1),'-')
        %pause        
        if nparts > mlvl
            tvec = tree_traverser(xcells,zcells,nparts,mlvl,ctgry,gcells,xbnds,zbnds,Ntrunc,gam,ep);            
        else
            tvec = near_neighbor_comp(xloc,zloc,xlist,zlist,gloc,glist,Ntrunc,gam,ep);            
        end
        Kvec(locinds,:) = Kfar + tvec;
    end
end