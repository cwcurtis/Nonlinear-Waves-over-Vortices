function Bendixson(xpos,zpos,gvals,n_bdry,Nvorts,markersize)
    
    hold on
    
    mc = max(abs(gvals));
    for kk = 1:Nvorts-n_bdry
        if(gvals(kk) ~= 0)
            a = min(gvals(kk)/mc,0);
            b = max(gvals(kk)/mc,0);
            plot(xpos(kk),zpos(kk),'.','MarkerSize',markersize','color',[1+a 1+a-b 1-b])
        end
    end
    
    plot(xpos(Nvorts-n_bdry+1:Nvorts),zpos(Nvorts-n_bdry+1:Nvorts),'k.','MarkerSize',1);
    hold off
    
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    xlim([-1 1])
    ylim([-1 2])
end