function Kloc = near_neighbor_comp_periodic(xloc,zloc,xfar,zfar,gnear,gfar,rval,Mx)

% and now the in-cell interactions
[xf,xl] = meshgrid(xloc,xloc);
[zf,zl] = meshgrid(zloc,zloc);
p2M = pi/(2*Mx);

dx = xl-xf;
dz = zl-zf;
dzp = zl+zf;
dzz = dx + 1i*dz + eye(size(dx));
dzzp = dx + 1i*dzp;

diff = abs(dzz).^2;
diffp = abs(dzzp).^2;

er = exp(-diff/(2*rval^2));
scv = (1-2*er).*er./(p2M*dzz);
erp = exp(-diffp/(2*rval^2));
scvp = (1-2*erp).*erp./(p2M*dzzp);

%tval = cot(p2M*dzz);
pker = cot(p2M*dzz) + scv;
pker = pker - diag(diag(pker));
pkerp = cot(p2M*dzzp) + scvp;
ptot = pker - pkerp;
Kvals = ptot*gnear;
Kloc1 = imag(Kvals);
Kloc2 = real(Kvals);

if ~isempty(xfar)
    [xf,xl] = meshgrid(xfar,xloc);
    [zf,zl] = meshgrid(zfar,zloc);
    dx = xl-xf;
    dz = zl-zf;
    dzp = zl+zf;
    dzz = dx + 1i*dz;
    dzzp = dx + 1i*dzp;

    diff = abs(dzz).^2;
    diffp = abs(dzzp).^2;

    er = exp(-diff/(2*rval^2));
    scv = (1-2*er).*er./(p2M*(dx+1i*dz));
    erp = exp(-diffp/(2*rval^2));
    scvp = (1-2*erp).*erp./(p2M*(dx+1i*dzp));

    pker = cot(p2M*dzz) + scv;
    pkerp = cot(p2M*dzzp) + scvp;
    ptot = pker - pkerp;
    Kvals = ptot*gfar;
    Kloc1 = Kloc1 + imag(Kvals);
    Kloc2 = Kloc2 + real(Kvals);
end


Kloc = [Kloc1 Kloc2];