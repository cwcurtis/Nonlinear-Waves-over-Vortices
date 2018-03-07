function Kfar = far_panel_exact_comp(xloc,zloc,xfar,zfar,gfar,gam)

nparts = length(xloc);
Kfar = zeros(nparts,2);

for jj=1:nparts
    
    dx = xloc(jj) - xfar;
    dzm = gam*(zloc(jj) - zfar);
    
    diff = dx.^2 + dzm.^2;
    kerx = -dzm./diff;
    kerz = dx./diff;   
    
    Kfar(jj,1) = sum(gfar.*kerx);
    Kfar(jj,2) = sum(gfar.*kerz);   
    
end    