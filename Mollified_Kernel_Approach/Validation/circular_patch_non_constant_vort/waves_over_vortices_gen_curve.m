function waves_over_vortices_gen_curve(Nx,mu,gam,omega,tf)
    close all
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    Rv = 1/25;
    zoff = .25;
    
    F = pi*omega*Rv^2/(4*gam);
    pmesh = linspace(0,2*pi,2*Nx);
     
    cp = cos(pmesh);
    sp = sin(pmesh);    
    
    disp('Froude number is:')
    disp(F)

    [xpos,zpos,gvals,ep,Nvorts] = initializer(Nx,omega,gam,Rv);
    simul_plot = 1; % Plot during computation.       0 - off, 1 - on
    n_bdry = 0;     % Number of points in cicular boundary.
    markersize = 10;
    
    dt = .1;
    nmax = round(tf/dt);
    
    xtrack = xpos;
    ztrack = zoff + zpos/gam;
    gtrack = gvals;
    
    inter = 2;
    plot_count = 1;
    no_of_evals = round(nmax/inter);
    times = zeros(no_of_evals,1);
    errors = zeros(no_of_evals,1);
    Vcnt = zeros(no_of_evals+1,1);
    Vcnt(1) = Nvorts;
    
    xellip = Rv*cp;
    yellip = Rv*sp/gam+zoff;
       
    u = [xpos;zpos]; %velocity vector field
    % Make folder
    S = make_folder(Nx/2,Nx,mu,gam,F,tf);
    filename = strcat(S, '/', '/waves_over_vortices.gif');
        
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
    clf
    if simul_plot
        figure(1)
        Bendixson(xpos,zoff + zpos/gam,gvals,n_bdry,Nvorts,markersize);
        hold on
            scatter(xellip,yellip,.1)
        axis equal
        hold off
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','DelayTime',0,'loopcount',inf);
    end
   
    tic
    for jj=1:nmax
      
        % Now update the vortex positions                                   
        u = vort_update_on_molly_non_periodic(mu,ep,u,gvals,Nvorts,dt);
            
        if(mod(jj,inter)==0)
            times(plot_count) = (jj-1)*dt;
            xpos = u(1:Nvorts);
            zpos = zoff + u(Nvorts+1:2*Nvorts)/gam;
            
            error = interp_error(Nx,xpos,zpos,gvals,ep,gam,omega,Rv,zoff,xellip,yellip);
            errors(plot_count) = error;
                        
            xtrack = [xtrack;xpos];
            ztrack = [ztrack;zpos];
            gtrack = [gtrack;gvals];
            plot_count = plot_count + 1;
            
            Vcnt(plot_count) = Nvorts;   
            
            if simul_plot
                figure(1)
                Bendixson(xpos,zpos,gvals,n_bdry,Nvorts,markersize);
                hold on
                    scatter(xellip,yellip,.1)
                    axis equal
                hold off        
                frame = getframe(1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                imwrite(imind,cm,filename,'gif','DelayTime',0,'writemode','append');
            end          
                        
            [xpud,zpud,gvud] = recircer(gvals,xpos,gam*(zpos-zoff),Nx);
            Nvorts = length(gvud);
            gvals = gvud;
            u = zeros(2*Nvorts,1);
            u(1:Nvorts) = xpud;
            u(Nvorts+1:2*Nvorts) = zpud;            
                        
        end       
    end
    
    toc
    
    figure(2)
    plot(times,errors,'k-','LineWidth',2)
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Animate waves over vortices
    if ~simul_plot
        figure(1)
        gif_my_gif(xtrack,ztrack,gtrack,n_bdry,Vcnt,plot_count,filename,markersize);
    end
    
end