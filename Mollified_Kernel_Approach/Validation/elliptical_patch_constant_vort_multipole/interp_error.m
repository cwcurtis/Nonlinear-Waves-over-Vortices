function error = interp_error(xpos,zpos,gvals,rval,gam,omega,av,bv,zoff)
    
    chi = @(r) 1./pi.*(2.*exp(-(r.^2))-.5*exp(-(r/sqrt(2)).^2));
    diffmat = sqrt((xpos*ones(1,length(xpos))-ones(length(xpos),1)*xpos').^2 + gam.^2*(zpos*ones(1,length(zpos))-ones(length(zpos),1)*zpos').^2)/rval;
    imat = gam*chi(diffmat)*gvals/(rval^2);
    %scatter3(xpos,zpos,imat)
    %pause
    rmat = omega*(1-(xpos/av).^2-gam^2*((zpos-zoff)/bv).^2).^3;
    inds = rmat < 0;
    rmat(inds) = 0;
    
    l2r = sqrt(sum(rmat.^2));
    error = sqrt(sum((rmat-imat).^2))/l2r;
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