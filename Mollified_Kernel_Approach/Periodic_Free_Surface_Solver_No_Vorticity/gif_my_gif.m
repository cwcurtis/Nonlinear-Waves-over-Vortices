function gif_my_gif(Xmesh,Mx,eta_plot,t,S)
    filename = strcat(S, '/', '/waves_over_vortices.gif');
    
    for ii = 1:t
        
        Bendixson(Xmesh,Mx,eta_plot(:,ii))
        
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        if ii == 1
            imwrite(imind,cm,filename,'gif','DelayTime',0,'loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',0,'writemode','append');
        end
    end
    
    imwrite(imind,cm,strcat(S,'.png'),'png');
end