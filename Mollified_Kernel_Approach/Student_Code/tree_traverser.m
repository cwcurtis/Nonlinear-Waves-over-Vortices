function Kvec = tree_traverser(xpos,zpos,nparts,mlvl,ctgry,xbnds,zbnds,gvals,Ntrunc,gam,ep)

Kvec = zeros(nparts,2);
ntot = length(xpos);

xmin = xbnds(1);
xmax = xbnds(2);
zmin = zbnds(1);
zmax = zbnds(2);

if strcmp(ctgry,'corner')
    ccnt = 4;
    rcnt = 4;
    nblcks = 16;
    jvals = [1;2;5;6];
elseif strcmp(ctgry,'hedge')
    ccnt = 6;
    rcnt = 4;
    nblcks = 24;
    jvals = [3;4;9;10];
elseif strcmp(ctgry,'vedge')
    ccnt = 4;
    rcnt = 6;
    nblcks = 24;
    jvals = [9;10;13;14];
elseif strcmp(ctgry,'interior')
    ccnt = 6;
    rcnt = 6;
    nblcks = 36;
    jvals = [15;16;21;22];
end

dx = (xmax-xmin)/ccnt;
dz = (zmax-zmin)/rcnt;
vvecs = zeros(ntot,ccnt);
hvecs = zeros(ntot,rcnt);

bvecs = zeros(ntot,nblcks);
ilists = zeros(ntot,nblcks);
farfield = zeros(ntot,nblcks);
   
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

cvals = zeros(3,1);

for ll=1:4
        
    jj = jvals(ll);
        
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
        
    npartsloc = sum(bvecs(:,jj));
    cvals(ll) = npartsloc;
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
    
    if npartsloc > 0
        Kfar = far_panel_comp(xloc,zloc,xfar,zfar,gfar,gam);        
        if npartsloc > mlvl
            if(ll==1) 
                Kvec(1:npartsloc,:) = Kfar + tree_traverser(xlist,zlist,npartsloc,mlvl,ctgry,xbnds,zbnds,gloc,Ntrunc,gam,ep);
            else
                Kvec(cvals(ll-1)+1:cvals(ll-1)+npartsloc,:) = Kfar + tree_traverser(xlist,zlist,npartsloc,mlvl,ctgry,xbnds,zbnds,gloc,Ntrunc,gam,ep); 
            end
        else
            cominds = (ilists(:,jj)-bvecs(:,jj))~=0;
            xfar = xpos(cominds);
            zfar = zpos(cominds);
            gnear = gvals(locinds);
            gfar = gvals(cominds);
            Kloc = near_neighbor_comp(xloc,zloc,xfar,zfar,gnear,gfar,Ntrunc,gam,ep);
            if(ll==1) 
                Kvec(1:npartsloc,:) = Kfar + Kloc;
            else
                Kvec(cvals(ll-1)+1:cvals(ll-1)+npartsloc,:) = Kfar + Kloc; 
            end
        end
    end
end