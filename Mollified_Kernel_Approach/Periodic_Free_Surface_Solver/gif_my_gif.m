function gif_my_gif(Xmesh,eta_plot,xtrack,ztrack,gtrack,n_bdry,Vcnt,t,filename,markersize)
    shft_ind = 0;
    for ii = 1:t
        
        Nvortl = Vcnt(ii);
        inds = (1+shft_ind):(Nvortl+shft_ind);
        Bendixson(Xmesh,eta_plot(:,ii),xtrack(inds),ztrack(inds),gtrack(inds),n_bdry,Nvortl,markersize)
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
end