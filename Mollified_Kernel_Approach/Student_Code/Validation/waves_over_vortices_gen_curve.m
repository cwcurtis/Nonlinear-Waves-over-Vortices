function waves_over_vortices_gen_curve(Nx,M,mu,gam,omega,ep,tf,Ntrunc)
    close all
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    Rv = 1/25;
    zoff = .25;
    axs = 2;
    cfun = @(x,z,axs,gam,Rv,zoff) 1/Rv^2*(x.^2/axs.^2 + (z-zoff).^2*gam^2);
    F = pi*omega*axs*Rv^2/gam;
    pmesh = linspace(0,2*pi,1e2);
    cp = cos(pmesh);
    sp = sin(pmesh);
    if axs>1
        avel = mu/gam*omega*axs/(axs+1)^2;
    else
        avel = 0;
    end
    disp('Effective Aspect Ratio is')
    disp(max([axs*gam 1/(axs*gam)]))
    [xpos,zpos,gvals,Gamma,Nvorts] = initializer(Nx,F,M,gam,cfun,Rv,zoff,axs);
    simul_plot = 0; % Plot during computation.       0 - off, 1 - on
    n_bdry = 0;     % Number of points in cicular boundary.
    markersize = 10;
    
    dt = 5e-2;
    nmax = round(tf/dt);
    
    xtrack = xpos;
    ztrack = zpos;
    gtrack = gvals;
    
    inter = 2;
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
    clf
    if simul_plot
        figure(1)
        Bendixson(xpos,zpos,gvals,n_bdry,Nvorts,markersize);
        xellip = axs*Rv*cp;
        yellip = Rv*sp/gam+zoff;
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
        u = vort_update_on_molly_non_periodic(mu,gam,ep,u,gvals,Nvorts,dt);
        
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
                ca = cos(avel*(jj-1)*dt);
                sa = sin(avel*(jj-1)*dt);
                xellip = Rv*(ca*axs*cp-sa*sp/gam);
                yellip = Rv*(ca*sp/gam + sa*axs*cp) + zoff;        
                hold on
                scatter(xellip,yellip,.1)
                axis equal
                hold off        
                frame = getframe(1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                imwrite(imind,cm,filename,'gif','DelayTime',0,'writemode','append');
            end          
            
            [xpud,zpud,gvud] = recircer(gvals,xpos,zpos,Nx);
            Nvorts = length(gvud);
            gvals = gvud;
            u = zeros(2*Nvorts,1);
            u(1:Nvorts) = xpud;
            u(Nvorts+1:2*Nvorts) = zpud;            
        end       
    end
    
    toc
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Animate waves over vortices
    if ~simul_plot
        figure(1)
        gif_my_gif(xtrack,ztrack,gtrack,n_bdry,Vcnt,avel,plot_count,filename,markersize);
    end
    
end