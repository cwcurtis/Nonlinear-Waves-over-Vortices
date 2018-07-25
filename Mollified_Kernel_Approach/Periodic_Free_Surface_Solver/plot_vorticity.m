function plot_vorticity(S,ep,gam,xpos,zpos,gvals)

    chi = @(r) 1./pi.*(2.*exp(-(r.^2))-.5*exp(-(r/sqrt(2)).^2));
    diffmat = sqrt((xpos*ones(1,length(xpos))-ones(length(xpos),1)*xpos').^2 + gam^2*(zpos*ones(1,length(zpos))-ones(length(zpos),1)*zpos').^2);
    ivec = gam*chi(diffmat/ep)*gvals/(ep^2);
    Finterp = scatteredInterpolant(xpos,zpos,ivec);
    xmin = min(xpos);
    xmax = max(xpos);
    ddx = (xmax-xmin)/100;
    zmin = min(zpos);
    zmax = max(zpos);
    ddz = (zmax-zmin)/100;
    [Xmm,Zmm] = meshgrid((xmin:ddx:xmax),(zmin:ddz:zmax));
    
    surface(Xmm,Zmm,Finterp(Xmm,Zmm),'LineStyle','none')
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$x$','Interpreter','LaTeX','FontSize',30)
    ylabel('$z$','Interpreter','LaTeX','FontSize',30)    
    
    savefig(strcat(S, '/', 'vorticity'))