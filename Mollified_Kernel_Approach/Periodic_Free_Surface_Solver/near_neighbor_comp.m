function Kloc = near_neighbor_comp(xloc,zloc,xfar,zfar,gnear,gfar,rval)

% and now the in-cell interactions
[xf,xl] = meshgrid(xloc,xloc);
[zf,zl] = meshgrid(zloc,zloc);
dx = xl-xf;
dz = zl-zf;
dzp = zl+zf;
dx2 = dx.^2;
diff = dx2 + dz.^2 + eye(size(dx));
diffp = dx2 + dzp.^2;
er = exp(-diff/(2*rval^2));
scv = (1-er).*(1+2*er)./diff;
erp = exp(-diffp/(2*rval^2));
scvp = (1-erp).*(1+2*erp)./diffp;

Kloc1 = (-dz.*scv + dzp.*scvp)*gnear;
Kloc2 = (dx.*(scv-scvp))*gnear;

if ~isempty(xfar)
    [xf,xl] = meshgrid(xfar,xloc);
    [zf,zl] = meshgrid(zfar,zloc);
    dx = xl-xf;
    dz = zl-zf;
    dzp = zl+zf;
    dx2 = dx.^2;
    diff = dx2 + dz.^2;
    diffp = dx2 + dzp.^2;
    er = exp(-diff/(2*rval^2));
    scv = (1-er).*(1+2*er)./diff;
    erp = exp(-diffp/(2*rval^2));
    scvp = (1-erp).*(1+2*erp)./diffp;
    Kloc1 = Kloc1 + (-dz.*scv + dzp.*scvp)*gfar;
    Kloc2 = Kloc2 + (dx.*(scv-scvp))*gfar;
end


Kloc = [Kloc1 Kloc2];