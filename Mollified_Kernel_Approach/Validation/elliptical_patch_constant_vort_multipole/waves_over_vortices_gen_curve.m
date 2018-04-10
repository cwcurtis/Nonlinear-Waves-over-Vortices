function waves_over_vortices_gen_curve(Nx,mu,gam,omega,tf)
    close all
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    Rv = 1/25;
    zoff = .35;
    axs = 2;
    av = Rv*axs;
    bv = Rv;
    
    F = pi*omega*av*bv/gam;
     
    [xpos,zpos,gvals,rval,Nvorts] = initializer(Nx,omega,gam,av,bv);
    
    dt = .1;
    nmax = round(tf/dt);
    
    zpos = zoff + zpos/gam;
    
    inter = 10;
    plot_count = 1;
    no_of_evals = round(nmax/inter);
    Vcnt = zeros(no_of_evals+1,1);
    Vcnt(1) = Nvorts;
    xtrack = xpos;
    ztrack = zpos;
    gtrack = gvals;
    
    u = [xpos;zpos]; %velocity vector field
    
    % Make folder
    S = make_folder(Nx/2,Nx,mu,gam,F,tf,av/bv);
    clf
    n_bdry = 0;     % Number of points in cicular boundary.
    markersize = 10;
        
    for jj=1:nmax      
        % Now update the vortex positions         
        u = vort_update_multipole(mu,gam,rval,u,gvals,Nvorts,dt);           
            
        if(mod(jj,inter)==0)
            plot_count = plot_count + 1;
            
            xtrack = [xtrack;xpos];
            ztrack = [ztrack;zpos];
            gtrack = [gtrack;gvals];
            Vcnt(plot_count) = Nvorts;            
        
            [xpos,zpos,gvals] = recircer(gvals,u(1:Nvorts),u(Nvorts+1:2*Nvorts),Nx);
            Nvorts = length(gvals);
            disp('Number of Vortices is')
            disp(Nvorts)
            u = zeros(2*Nvorts,1);
            u(1:Nvorts) = xpos;
            u(Nvorts+1:2*Nvorts) = zpos;                    
        end
    end     
     
    figure(1)
    gif_my_gif(xtrack,ztrack,gtrack,n_bdry,Vcnt,plot_count,S,markersize);    
end