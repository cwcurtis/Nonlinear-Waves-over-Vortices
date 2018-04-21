function Bendixson(Xmesh,Mx,eta)
    plot(Xmesh,eta+1,'k','LineWidth',2)
    
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('x','Interpreter','LaTeX','FontSize',30)
    ylabel('z','Interpreter','LaTeX','FontSize',30)
    xlim([-Mx Mx])
    ylim([0 1.5])
end