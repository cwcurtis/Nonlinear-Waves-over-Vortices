function Kloc = near_neighbor_comp(xloc,zloc,xfar,zfar,gnear,gfar,gam,rval)

nparts = length(xloc);
Kloc = zeros(nparts,2);

for jj=1:nparts
    
    dx = xloc(jj) - xfar;
    dzm = gam*(zloc(jj) - zfar);
    
    diff = dx.^2 + dzm.^2;
    er = exp(-diff/(2*rval^2));
    scv = (1-er).*(1+2*er)./diff;
    kerx = -dzm.*scv;
    kerz = dx.*scv;   
    
    Kloc(jj,1) = sum(gfar.*kerx);
    Kloc(jj,2) = sum(gfar.*kerz);     
    
end    

for jj=1:nparts
    
    dx = xloc(jj) - xloc;
    dzm = gam*(zloc(jj) - zloc);
    
    diff = dx.^2 + dzm.^2;
    er = exp(-diff/(2*rval^2));
    scv = (1-er).*(1+2*er)./diff;
    kerx = -dzm.*scv;
    kerz = dx.*scv;   
    
    Kloc(jj,1) = Kloc(jj,1) + sum(gnear.*kerx);
    Kloc(jj,2) = Kloc(jj,2) + sum(gnear.*kerz);     
    
end    