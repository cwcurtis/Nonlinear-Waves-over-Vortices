function qvals = far_panel_comp(xfar,zfar,gfar,gam,xc,zc,pval)
c = xc + 1i*gam*zc;
zcf = c - (xfar+1i*gam*zfar);
qvals = zeros(pval+2,1);
qvals(1) = sum(gfar);
zi = gfar;
for jj=1:pval+1   
    %qvals(jj+1) = (-1)^(jj-1)*sum(zcf.^(jj-1).*gfar);
    qvals(jj+1) = sum(zi);
    zi = -zcf.*zi;
end



