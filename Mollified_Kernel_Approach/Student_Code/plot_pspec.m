function plot_pspec(K,S,pspec,eta0)
%    plot(-K+1:K,pspec,'k',-K+1:K,eta0,'LineWidth',2)
    plot(-K:K-1,pspec,'k','LineWidth',2)
    
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('k','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylabel('\eta(k,t)','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylim([-20 0])
    xlim([-K K-1])
    
    savefig(strcat(S, '/', 'fig2'))
end