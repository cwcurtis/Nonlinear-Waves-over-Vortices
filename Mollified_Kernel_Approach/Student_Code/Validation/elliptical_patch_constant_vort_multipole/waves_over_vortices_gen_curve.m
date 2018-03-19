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
    pmesh = linspace(0,2*pi,2*Nx);
     
    cp = cos(pmesh);
    sp = sin(pmesh);    
    
    [xpos,zpos,gvals,rval,Nvorts] = initializer(Nx,omega,gam,av,bv);
    
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
    
    u0 = [xpos;zpos]; %velocity vector field
        
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
    
     for jj=1:nmax
      
        % Now update the vortex positions         
        tic
        ut = vort_update_on_molly_non_periodic(mu,1,rval,u0,gvals,Nvorts,dt);
        toc
        
        tic
        u = vort_update_multipole(mu,1,rval,u0,gvals,Nvorts,dt);   
        toc
        disp(norm(u(1:Nvorts)-ut(1:Nvorts))/norm(ut(1:Nvorts)))
        disp(norm(u(Nvorts+1:2*Nvorts)-ut(Nvorts+1:2*Nvorts))/norm(ut(Nvorts+1:2*Nvorts)))
        
        %clf 
        
        %figure(1)
        %scatter(u(1:Nvorts),u(Nvorts+1:2*Nvorts),markersize,'MarkerFaceColor','k')
        
        %figure(2)
        %scatter(ut(1:Nvorts),ut(Nvorts+1:2*Nvorts),markersize,'MarkerFaceColor','r')
        
        %pause
        
        if(mod(jj,inter)==0)
            times(plot_count) = (jj-1)*dt;
            xpos = u(1:Nvorts);
            zpos = zoff + u(Nvorts+1:2*Nvorts)/gam;
            
            tv = mu/gam*omega*av*bv/(av+bv)^2*(jj-1)*dt;
            ca = cos(tv);
            sa = sin(tv);
            xxs = cp*av;
            yxs = sp*bv;
            xellip = xxs*ca - yxs*sa;
            yellip = (xxs*sa + yxs*ca)/gam + zoff;
            
            error = interp_error(xpos,zpos,gvals,rval,gam,omega,av,bv,tv,zoff);
            errors(plot_count) = error;
                        
            xtrack = [xtrack;xpos];
            ztrack = [ztrack;zpos];
            gtrack = [gtrack;gvals];
            plot_count = plot_count + 1;
            
            Vcnt(plot_count) = Nvorts;   
        
        end
        
        if(mod(jj,2*inter)==0)
            [xpud,zpud,gvud] = recircer(gvals,xpos,gam*(zpos-zoff),Nx);
            Nvorts = length(gvud);
            gvals = gvud;
            u = zeros(2*Nvorts,1);
            u(1:Nvorts) = xpud;
            u(Nvorts+1:2*Nvorts) = zpud;            
        end
            
    end
    
    %toc
    
    %figure(2)
    %plot(times,errors,'k-','LineWidth',2)
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
    
end