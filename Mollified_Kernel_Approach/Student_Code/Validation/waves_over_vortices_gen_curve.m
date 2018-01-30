function waves_over_vortices_gen_curve(Nx,mu,gam,omega,ep,tf,Ntrunc)
    close all
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    Rv = 1;
    zoff = .5;
    cfun = @(x,z) (1/Rv^2)*(x.^2 + gam^2*(z-zoff).^2);
    F = pi*Rv^8/(4*mu)*omega;
    disp('Froude number is')
    disp(F)
    [xpos,zpos,gvals,Gamma] = initializer(Nx,omega,cfun,gam,Rv,zoff);
    Nvorts = length(gvals);
    disp('Number of Starting Vortices is:')
    disp(Nvorts)
    simul_plot = 0; % Plot during computation.       0 - off, 1 - on
    n_bdry = 0;     % Number of points in cicular boundary.
    markersize = 10;
    
    dt = .05;
    nmax = round(tf/dt);
    
    xtrack = xpos;
    ztrack = zpos;
    gtrack = gvals;
    
    inter = 10;
    plot_count = 1;
    no_of_evals = round(nmax/inter);
    times = zeros(no_of_evals,1);
    Vcnt = zeros(no_of_evals+1,1);
    Vcnt(1) = Nvorts;
    
    u = [xpos;zpos]; %velocity vector field
    
    % Make folder
    S = make_folder(Nx/2,Nx,mu,gam,F,tf,Gamma);
    filename = strcat(S, '/', '/waves_over_vortices.gif');
        
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
    
    if simul_plot
        figure(1)
        Bendixson(xpos,zpos,gvals,n_bdry,Nvorts,markersize);
        
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','DelayTime',0,'loopcount',inf);
    end
   
    tic
    for jj=1:nmax
      
        % Now update the vortex positions                                   
        u = vort_update_on_molly_fourier(gam,mu,ep,F,u,gvals,Nvorts,Ntrunc,dt);
        
        if(mod(jj,inter)==0)
            times(plot_count) = (jj-1)*dt;
            plot_count = plot_count + 1;
            xpos = u(1:Nvorts);
            zpos = u(Nvorts+1:2*Nvorts);
            
            xtrack = [xtrack;xpos];
            ztrack = [ztrack;zpos];
            gtrack = [gtrack;gvals];
            Vcnt(plot_count) = Nvorts;   
            
            if simul_plot
                figure(1)
                Bendixson(xpos,zpos,gvals,n_bdry,Nvorts,markersize);
                frame = getframe(1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                imwrite(imind,cm,filename,'gif','DelayTime',0,'writemode','append');
            end          
            
        end
        
        if(mod(jj,10*inter)==0)        
            xpos = u(1:Nvorts);
            zpos = u(Nvorts+1:2*Nvorts);
            
            [xpos,zpos,gvals] = recircer(gvals,xpos,zpos,Nx);
            Nvorts = length(xpos);
            disp('Current number of vortices is')
            disp(Nvorts)
            u(1:Nvorts) = xpos;
            u(Nvorts+1:2*Nvorts) = zpos;
        end       
    end
    
    toc
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Animate waves over vortices
    if ~simul_plot
        figure(1)
        gif_my_gif(xtrack,ztrack,gtrack,n_bdry,Vcnt,plot_count,filename,markersize);
    end
    
end