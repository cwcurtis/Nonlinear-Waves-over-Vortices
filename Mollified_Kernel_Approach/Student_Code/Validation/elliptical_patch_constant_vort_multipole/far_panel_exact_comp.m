function Kfar = far_panel_exact_comp(xloc,zloc,xfar,zfar,gfar,rval)

nparts = length(xloc);
Kfar = zeros(nparts,2);

[xf,xl] = meshgrid(xfar,xloc);
[zf,zl] = meshgrid(zfar,zloc);
dx = xl-xf;
dz = zl-zf;
%dx = bsxfun(@minus,xloc,xfar.');
%dz = bsxfun(@minus,zloc,zfar.');
diff = dx.^2 + dz.^2;
er = exp(-diff/(2*rval^2));
scv = (1-er).*(1+2*er)./diff;
Kfar(:,1) = (-dz.*scv)*gfar;
Kfar(:,2) = (dx.*scv)*gfar;
