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
     
    [xpos,zpos,gvals,rval,Nvorts] = initializer(Nx,omega,gam,av,bv);
    
    dt = .1;
    nmax = round(tf/dt);
    
    zpos = zoff + zpos/gam;
    
    inter = 10;
    Nvortst = Nvorts;
    gvalst = gvals;
    ut = [xpos;zpos]; %velocity vector field
    u = [xpos;zpos]; %velocity vector field
    
     for jj=1:nmax      
        % Now update the vortex positions         
        
        tic
        u = vort_update_multipole(mu,gam,rval,u,gvals,Nvorts,dt);           
        toc
        
        tic
        ut = vort_update_on_molly_non_periodic(mu,gam,rval,ut,gvalst,Nvortst,dt);           
        toc
        
        disp(norm(ut-u)/norm(ut))
        
        if(mod(jj,inter)==0)
            [xpost,zpost,gvalst] = recircer(gvalst,ut(1:Nvortst),ut(Nvortst+1:2*Nvortst),Nx);
            Nvortst = length(gvalst);
            disp('Number of Vortices is')
            disp(Nvortst)
            ut = zeros(2*Nvortst,1);
            ut(1:Nvortst) = xpost;
            ut(Nvortst+1:2*Nvortst) = zpost;            
        end
                
        if(mod(jj,inter)==0)
            [xpos,zpos,gvals] = recircer(gvals,u(1:Nvorts),u(Nvorts+1:2*Nvorts),Nx);
            Nvorts = length(gvals);
            disp('Number of Vortices is')
            disp(Nvorts)
            u = zeros(2*Nvorts,1);
            u(1:Nvorts) = xpos;
            u(Nvorts+1:2*Nvorts) = zpos;                    
        end
     end     
end