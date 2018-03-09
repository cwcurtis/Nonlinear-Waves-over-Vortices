function qvals = far_panel_comp(xfar,zfar,gfar,gam,xc,zc,pval)

zcf = (xc-xfar) + 1i*gam*(zc-zfar);
nparts = length(zcf);
zmat = ones(nparts,pval+1);

for jj=2:pval+1   
    zmat(:,jj) = zmat(:,jj-1).*(-zcf);    
end

qvals = zmat.'*gfar;

