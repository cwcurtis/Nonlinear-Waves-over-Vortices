function Kvec = multi_pole_kernel(xpos,zpos,gvals,Ntrunc,gam,ep)

% Build kd-tree structure from xpos,zpos
Nvorts = length(xpos);
mlvl = floor(log2(Nvorts));
Kvec = zeros(Nvorts,2);

xmin = min(xpos);
xmax = max(xpos);
zmin = min(zpos);
zmax = max(zpos);

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
    xr = xmin + jj*dx;
    vvecs(:,jj) = (xpos>=xl).*(xpos<xr);    
end

for jj=1:rcnt
    zt = zmax - (jj-1)*dz;
    zb = zmax - jj*dz;
    hvecs(:,jj) = (zpos>=zb).*(zpos<zt); 
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
    
    xl = xmin + (col-lshift+1)*dx;
    xr = xmin + (col+rshift+1)*dx;
    zt = zmax - (row-dshift+1)*dz;
    zb = zmax - (row+ushift+1)*dz;
    xbnds = [xl;xr];
    zbnds = [zb;zt];
    
    nparts = sum(bvecs(:,jj));
    farfield(:,jj) = 1 - ilists(:,jj);
    
    locinds = bvecs(:,jj)~=0;
    lisinds = ilists(:,jj)~=0;
    farinds = farfield(:,jj)~=0;
    
    xloc = xpos(locinds);
    xlist = xpos(lisinds);
    xfar = xpos(farinds);
    
    zloc = zpos(locinds);
    zlist = zpos(lisinds);
    zfar = zpos(farinds);
    
    gloc = gvals(lisinds);
    gfar = gvals(farinds);
    
    if nparts > 0
        Kfar = far_panel_comp(xloc,zloc,xfar,zfar,gfar,gam);
        if nparts > mlvl
            Kvec(locinds,:) = Kfar + tree_traverser(xlist,zlist,nparts,mlvl,ctgry,xbnds,zbnds,gloc,Ntrunc,gam,ep);        
        else
            cominds = (ilists(:,jj)-bvecs(:,jj))~=0;
            xfar = xpos(cominds);
            zfar = zpos(cominds);
            gnear = gvals(locinds);
            gfar = gvals(cominds);
            Kloc = near_neighbor_comp(xloc,zloc,xfar,zfar,gnear,gfar,Ntrunc,gam,ep);
            Kvec(locinds,:) = Kfar + Kloc;
        end
    end
end