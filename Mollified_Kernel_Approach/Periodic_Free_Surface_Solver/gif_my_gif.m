function gif_my_gif(Xmesh,Mx,eta_plot,etanv_plot,xtrack,ztrack,gtrack,Vcnt,t,S,markersize)
    shft_ind = 0;
    filename = strcat(S, '/', '/waves_over_vortices.gif');
    for ii = 1:t
        
        Nvortl = Vcnt(ii);
        inds = (1+shft_ind):(Nvortl+shft_ind);
        Bendixson(Xmesh,Mx,eta_plot(:,ii),etanv_plot,xtrack(inds),ztrack(inds),gtrack(inds),markersize)
        shft_ind = shft_ind+Nvortl;
        
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