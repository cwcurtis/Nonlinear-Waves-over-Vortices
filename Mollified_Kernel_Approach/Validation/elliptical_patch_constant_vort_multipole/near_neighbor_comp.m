function Kloc = near_neighbor_comp(xloc,zloc,xfar,zfar,gnear,gfar,rval)

% and now the in-cell interactions
[xf,xl] = meshgrid(xloc,xloc);
[zf,zl] = meshgrid(zloc,zloc);
dx = xl-xf;
dz = zl-zf;
%dx = bsxfun(@minus,xloc,xloc');
%dz = bsxfun(@minus,zloc,zloc');
diff = dx.^2 + dz.^2 + eye(size(dx));
er = exp(-diff/(2*rval^2));
scv = (1-er).*(1+2*er)./diff;
Kloc1 = (-dz.*scv)*gnear;
Kloc2 = (dx.*scv)*gnear;

if ~isempty(xfar)
    [xf,xl] = meshgrid(xfar,xloc);
    [zf,zl] = meshgrid(zfar,zloc);
    dx = xl-xf;
    dz = zl-zf;
    %dx = bsxfun(@minus,xloc,xfar.');
    %dz = bsxfun(@minus,zloc,zfar.');
    diff = dx.^2 + dz.^2;
    er = exp(-diff/(2*rval^2));
    scv = (1-er).*(1+2*er)./diff;
    Kloc1 = Kloc1 + (-dz.*scv)*gfar;
    Kloc2 = Kloc2 + (dx.*scv)*gfar;
end


Kloc = [Kloc1 Kloc2];