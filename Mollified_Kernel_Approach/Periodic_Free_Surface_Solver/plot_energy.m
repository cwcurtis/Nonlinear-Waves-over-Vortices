function plot_energy(t,S,times,p1)
    %plot(times(1:t-1),p1(1:t-1),'k',times(1:t-1),p2(1:t-1),'k--',times(1:t-1),p3(1:t-1),'k-.','LineWidth',2)
    plot(times(1:t-1),p1(1:t-1),'k','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$t$','Interpreter','LaTeX','FontSize',30)
    ylabel('$\delta E(t)$','Interpreter','LaTeX','FontSize',30)
    
    %{
    % get ylim
    ym=min(abs(p1(1:t-1)));  
    
    % get order of magnitude
    e=log10(abs(ym));
    e=sign(e)*floor(abs(e));
    
    % get and rescale yticks
    yt=get(gca,'ytick')/10^e;
    % create tick labels
    ytl=cell(size(yt));
    for j=1:length(yt)
        % the space after the percent gives the same size to positive and
        % negative numbers. The number of decimal digits can be changed.
        ytl{j}=sprintf('% 3.2g',yt(j));
    end
    % set tick labels
    set(gca,'yticklabel',ytl);
    % place order of magnitude
    fs = get(gca,'fontsize');
    set(gca,'units','normalized');
    xl = xlim;
    yl = ylim;
    text(xl(1),yl(2),sprintf('x 10^{%d}',e),...
    'fontsize',fs,'VerticalAlignment','bottom');
    %}
    
    savefig(strcat(S, '/', 'energy'))
end