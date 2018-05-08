function plot_energy(t,S,times,p1,p2,p3)
    plot(times(1:t-1),p1(1:t-1),'k',times(1:t-1),p2(1:t-1),'k--',times(1:t-1),p3(1:t-1),'k-.','LineWidth',2)
    
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$t$','Interpreter','LaTeX','FontSize',30)
    savefig(strcat(S, '/', 'energy'))
end