function error = interp_error(Nx,xpos,zpos,gvals,ep,gam,om,Rv,zoff,xellip,yellip)

    xvals = linspace(-1,1,Nx)';
    zvals = linspace(0,1,Nx)';
    
    chi = @(r) 1/(ep^2*pi)*(2.*exp(-(r/ep).^2)-.5*exp(-(r/(ep*sqrt(2))).^2));
    omfunc = @(xj,yj) om*(1 - (xj.^2 + gam^2*(yj-zoff).^2)/Rv.^2).^3;
    
    imat = zeros(length(zvals),length(xvals));
    rmat = zeros(length(zvals),length(xvals));
    
    for jj=1:length(xvals)
       for kk=1:length(zvals)
          rval = sqrt((xvals(jj)-xpos).^2 + gam.^2*(zvals(kk)-zpos).^2);
          imat(kk,jj) = gam*sum(gvals.*chi(rval));
       end
    end
    
    ymax = max(yellip);
    ymin = min(yellip);
    jmx = 1 + floor(ymax*Nx);
    jmn = 1 + ceil(ymin*Nx);
    
    for ll = jmn:jmx
       ku = yellip < zvals(ll)+1/(2*Nx);
       kd = yellip >= zvals(ll)-1/(2*Nx);
       kvals = logical(ku.*kd);
       xl = min(xellip(kvals));
       xr = max(xellip(kvals));
       
       xvl = xvals >= xl;
       xvr = xvals <= xr;
       inds = logical(xvl.*xvr);
       rmat(ll,inds) = omfunc(xvals(inds),zvals(ll)); 
    end
    
    %{
    figure(3)
    surf(xvals,zvals,log10(abs(rmat-imat)),'LineStyle','none')
    
    figure(4)
    surf(xvals,zvals,imat,'LineStyle','none')
    
    figure(5)
    surf(xvals,zvals,rmat,'LineStyle','none')
        
    pause
    %}
    
    l2r = sqrt(2*(1/Nx).^2*sum(sum(rmat.^2)));
    error = sqrt(2*(1/Nx).^2*sum(sum((rmat-imat).^2)))/l2r;
    