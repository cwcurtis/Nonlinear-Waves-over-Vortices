function qvals = far_panel_comp(xfar,zfar,gfar,xc,zc,pval,Mx)

%zcf = (xc-xfar) + 1i*(zc-zfar);
zcf = tan(pi/(2*Mx)*((xc-xfar) + 1i*(zc-zfar)));
nparts = length(zcf);
%zmat = ones(nparts,pval+1);
zmat = ones(nparts,pval+2);

for jj=2:pval+2   
    zmat(:,jj) = zmat(:,jj-1).*zcf;        
end

qvals = zmat.'*gfar;

