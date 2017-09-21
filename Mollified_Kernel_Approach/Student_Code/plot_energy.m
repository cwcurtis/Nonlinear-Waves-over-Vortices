function plot_energy(t,S,times,p1,p2,p3)
    plot(times(1:t-1),p1(1:t-1),'k',times(1:t-1),p2(1:t-1),'k--',times(1:t-1),p3(1:t-1),'k-.','LineWidth',2)
    
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('t','FontName','Helvetica','FontSize',30,'FontWeight','bold')
%    ylabel('E(t)','FontName','Helvetica','FontSize',30,'FontWeight','bold')

    savefig(strcat(S, '/', 'fig3'))
end