function plot_com(S,comt,comx)

    plot(comt,comx,'k','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$t$','Interpreter','LaTeX','FontSize',30)
    ylabel('$x_{c}(t)$','Interpreter','LaTeX','FontSize',30)    
    
    savefig(strcat(S, '/', 'com'))