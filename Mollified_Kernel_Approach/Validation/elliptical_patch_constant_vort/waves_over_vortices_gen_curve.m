function waves_over_vortices_gen_curve(Nx,mu,gam,omega,tf)
    close all
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    Rv = 1/25;
    zoff = .25;
    axs = 2;
    av = Rv*axs;
    bv = Rv;
    
    F = pi*omega*av*bv/gam;
    Gamma = F;
    pmesh = linspace(0,2*pi,2*Nx);
     
    cp = cos(pmesh);
    sp = sin(pmesh);    
    
    [xpos,zpos,gvals,rval,Nvorts] = initializer(Nx,omega,gam,av,bv);
    simul_plot = 0; % Plot during computation.       0 - off, 1 - on
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
        Bendixson(xpos,zoff + zpos/gam,gvals,n_bdry,Nvorts,markersize);
        xellip = av*cp;
        yellip = bv*sp/gam+zoff;
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
        u = gpuArray(u);
        gvals = gpuArray(gvals);
    
        u = vort_update_on_molly_non_periodic(mu,rval,u,gvals,Nvorts,dt);
        
        u = gather(u);
        gvals = gather(gvals);
        
        if(mod(jj,inter)==0)
            times(plot_count) = (jj-1)*dt;
            xpos = u(1:Nvorts);
            zpos = zoff + u(Nvorts+1:2*Nvorts)/gam;
            
            %tv = mu/gam*omega*av*bv/(av+bv)^2*(jj-1)*dt;
            %ca = cos(tv);
            %sa = sin(tv);
            %xxs = cp*av;
            %yxs = sp*bv;
            %xellip = xxs*ca - yxs*sa;
            %yellip = (xxs*sa + yxs*ca)/gam + zoff;
            
            %error = interp_error(xpos,zpos,gvals,rval,gam,omega,av,bv,tv,zoff);
            %errors(plot_count) = error;
                        
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
            
        end 
                
        if(mod(jj,2*inter)==0)
            [xpud,zpud,gvud] = recircer(gvals,xpos,gam*(zpos-zoff),Nx);
            Nvorts = length(gvud);
            disp(Nvorts)
            gvals = gvud;
            u = zeros(2*Nvorts,1);
            u(1:Nvorts) = xpud;
            u(Nvorts+1:2*Nvorts) = zpud;            
        end
            
    end
    
    toc
    
    %figure(2)
    %plot(times,errors,'k-','LineWidth',2)
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Animate waves over vortices
    if ~simul_plot
        figure(1)
        gif_my_gif(xtrack,ztrack,gtrack,n_bdry,Vcnt,plot_count,filename,markersize);
    end
    
end