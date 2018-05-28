function plot_pspec(K,S,pspec,eta0)
    plot(-K+1:K,pspec,'k',-K+1:K,eta0,'b--','LineWidth',2)
    %plot(-K:K-1,pspec,'k','LineWidth',2)
    
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$k$','Interpreter','LaTeX','FontSize',30)
    ylabel('$\log_{10}\left|\hat{\eta}(k,t_{f})\right|$','Interpreter','LaTeX','FontSize',30)
    ylim([-20 0])
    xlim([-K K-1])
    
    savefig(strcat(S, '/', 'spectrum'))
end