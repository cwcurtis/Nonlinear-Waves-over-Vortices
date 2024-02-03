function Bendixson(xpos,zpos,gval,r0,n_bdry,Nvorts,markersize)
    clf
    hold on
    
    mc = max(abs(gval));
    for kk = 1:Nvorts-n_bdry
        if(gval ~= 0)
            a = min(gval/mc,0);
            b = max(gval/mc,0);
            plot(xpos(kk),zpos(kk),'.','MarkerSize',markersize','color',[1+a 1+a-b 1-b])
        end
    end
    
    plot(xpos(Nvorts-n_bdry+1:Nvorts),zpos(Nvorts-n_bdry+1:Nvorts),'k.','MarkerSize',1);
    
    hold off
    axis equal
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    xlim([-1.1*r0 1.1*r0])
    ylim([-1.1*r0 1.1*r0])
end