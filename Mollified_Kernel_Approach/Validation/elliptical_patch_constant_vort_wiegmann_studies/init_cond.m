function phix = init_cond(Xmesh,gam,xpos,zpos,gvals)

phix = zeros(length(Xmesh),1);

fphix = @(x,zj) sinh(pi*gam*zj).*(cosh(pi*gam*zj) - cosh(pi*gam).*cos(pi*x))./( (cosh(pi*gam*(1-zj)) - cos(pi*x)).*( cosh(pi*gam*(1+zj)) - cos(pi*x) ) );

for jj=1:length(Xmesh)

    xvec = 1/2*gvals.*( fphix(Xmesh(jj)-xpos,zpos) );
    
    phix(jj) = sum( xvec );    
    
end
