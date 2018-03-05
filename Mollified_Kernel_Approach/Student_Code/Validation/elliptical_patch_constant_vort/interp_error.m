function error = interp_error(xpos,zpos,gvals,rval,gam,om,av,bv,tv,zoff)


    Nx = 800;
    xvals = linspace(-.09,.09,Nx)';
    zvals = linspace(.12,.38,Nx)';
    
    ct = cos(tv);
    st = sin(tv);
    ax = ct^2/av^2 + st^2/bv^2;
    bx = 2*gam*(1/av^2 - 1/bv^2)*ct*st;
    cx = gam^2*(ct^2/bv^2 + st^2/av^2);
    
    chi = @(r) 1./pi.*(2.*exp(-(r.^2))-.5*exp(-(r/sqrt(2)).^2));
    
    imat = zeros(length(zvals),length(xvals));
    rmat = zeros(length(zvals),length(xvals));
    cfun = @(x,z) (x/av).^2 + (z/bv).^2;
    ofun = @(x,z) om*exp(log(eps)*(cfun(x,z)).^18);
    
    for jj=1:length(xvals)
       for kk=1:length(zvals)
          diff = sqrt((xvals(jj)-xpos).^2 + gam.^2*(zvals(kk)-zpos).^2)/rval;          
          imat(kk,jj) = gam*sum(gvals.*chi(diff)/(rval^2));
       end
    end
    
    for ll = 1:Nx
        inds = ax*xvals.^2 + bx*xvals*(zvals(ll)-zoff) + cx*(zvals(ll)-zoff)^2 <= 1;
        if length(inds)>1
            rmat(ll,inds) = om; 
        end
    end
    
    l2r = sqrt(sum(sum(rmat.^2)));
    error = sqrt(sum(sum((rmat-imat).^2)))/l2r;
        
    figure(3)
    surf(xvals,zvals,log10(abs(rmat-imat)),'LineStyle','none')
    
    figure(4)
    surf(xvals,zvals,imat,'LineStyle','none')
    
    figure(5)
    surf(xvals,zvals,rmat,'LineStyle','none')
     
    disp(error)
    
    pause
       
    %error = max(max(abs(rmat-imat)))/max(max(abs(rmat)));