function Bendixson(Xmesh,Mx,eta,eta_nv,xpos,zpos,gvals,markersize)
    
    clf
    hold on
    
    plot(Xmesh,eta+1,'k','LineWidth',2)
    mgv = max(abs(gvals));
    
    pgsinds = gvals > 0;
    gps = gvals(pgsinds)/mgv;
    xps = xpos(pgsinds);
    zps = zpos(pgsinds);    
    scatter(xps,zps,markersize,[ones(length(gps),1) 1-gps 1-gps],'filled')
    
    hold off
    
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('x','Interpreter','LaTeX','FontSize',30)
    ylabel('z','Interpreter','LaTeX','FontSize',30)
    xlim([-Mx Mx])
    ylim([0 1.5])
end