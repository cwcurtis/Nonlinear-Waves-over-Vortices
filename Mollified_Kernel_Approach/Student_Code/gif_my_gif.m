function gif_my_gif(Xmesh,eta_plot,xtrack,ztrack,gvals,n_bdry,Nvorts,inter,t,filename,markersize)
    for ii = 1:t
        Bendixson(Xmesh,eta_plot(:,ii),xtrack(:,floor((ii-1)*inter)+1),ztrack(:,floor((ii-1)*inter)+1),gvals,n_bdry,Nvorts,markersize)
        
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        if ii == 1
            imwrite(imind,cm,filename,'gif','DelayTime',0,'loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',0,'writemode','append');
        end
    end
end