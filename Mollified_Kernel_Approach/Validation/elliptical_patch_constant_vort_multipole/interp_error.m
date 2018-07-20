function error = interp_error(xpos,zpos,gvals,rval,gam,om,av,bv,zoff,rmat)

    Nx = 200;
    xvals = linspace(-.2,.2,Nx)';
    zvals = linspace(.1,.9,Nx)';   
    
    chi = @(r) 1./pi.*(2.*exp(-(r.^2))-.5*exp(-(r/sqrt(2)).^2));
    imat = zeros(length(zvals),length(xvals));
    
    for jj=1:length(xvals)
       for kk=1:length(zvals)
          diff = sqrt((xvals(jj)-xpos).^2 + gam.^2*(zvals(kk)-zpos).^2)/rval;          
          imat(kk,jj) = gam*sum(gvals.*chi(diff)/(rval^2));
       end
    end   
    
    l2r = sqrt(sum(sum(rmat.^2)));
    error = sqrt(sum(sum((rmat-imat).^2)))/l2r;
    %{
    disp(error)
        
    figure(3)
    surf(xvals,zvals,imat,'LineStyle','none')
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$x$','Interpreter','LaTeX','FontSize',30)
    ylabel('$z$','Interpreter','LaTeX','FontSize',30)
    
    figure(4)
    plot(xvals,imat(Nx/2+1,:),'LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$x$','Interpreter','LaTeX','FontSize',30)
    ylabel('$z$','Interpreter','LaTeX','FontSize',30)
    %}       