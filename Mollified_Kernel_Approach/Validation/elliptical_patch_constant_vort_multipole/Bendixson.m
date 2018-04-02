function Bendixson(xpos,zpos,gvals,markersize)
    clf 
    hold on
    mgv = max(abs(gvals));
    
    pgsinds = gvals > 0;
    gps = gvals(pgsinds)/mgv;
    xps = xpos(pgsinds);
    zps = zpos(pgsinds);    
    scatter(xps,zps,markersize,[ones(length(gps),1) 1-gps 1-gps],'filled')
    
    %{
    ngsinds = gvals < 0;
    gns = abs(gvals(ngsinds))/mgv;
    xns = xpos(ngsinds);
    zns = zpos(ngsinds);        
    scatter(xns,zns,markersize,[1-gns 1-gns ones(length(gns),1)],'filled')
    %}
    
    hold off
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('x','Interpreter','LaTeX','FontSize',30)
    ylabel('z','Interpreter','LaTeX','FontSize',30)
    xlim([-1 1])
    ylim([0 1])
end