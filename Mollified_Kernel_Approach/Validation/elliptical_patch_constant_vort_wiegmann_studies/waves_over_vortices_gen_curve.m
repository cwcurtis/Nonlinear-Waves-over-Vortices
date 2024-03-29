function waves_over_vortices_gen_curve(Nvorts,omega,gval,tf)
    close all
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    
    [xpos,zpos,r0,lscl] = initializer(omega,gval,Nvorts);
    
    Nvorts = size(xpos,2);
    disp(Nvorts)
    ndists = floor(sqrt(8*pi)*r0/lscl);

    n_bdry = 0;     % Number of points in cicular boundary.
    markersize = 10;
    
    dt = .1;
    nmax = round(tf/dt);
    
    inter = 10;
    plot_count = 1;
    no_of_evals = round(nmax/inter);
    times = zeros(no_of_evals,1);
    ang_moment = zeros(no_of_evals,1);
    hamiltonian = zeros(no_of_evals,1);
    radial_densities = zeros(no_of_evals, ndists);

    u = [xpos';zpos']; %velocity vector field
    
    % Make folder
    %S = make_folder(Nx,Nx,tf,Gamma);
    %filename = strcat(S, '/', '/waves_over_vortices.gif');
        
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
    
    tic
    
    for jj=1:nmax
      
        % Now update the vortex positions                                   
        %u = gpuArray(u);
        %gvals = gpuArray(gvals);
    
        u = vort_update_on_molly_non_periodic(u,omega,gval,Nvorts,dt);
        %u = dmm(u,dt,rval,omega,gvals,Nvorts);        
        
        %u = gather(u);
        %gvals = gather(gvals);
        
        if(mod(jj,inter)==0)
            times(plot_count) = (jj-1)*dt;
            xpos = u(1:Nvorts);
            zpos = u(Nvorts+1:2*Nvorts);
            
            ang_moment(plot_count) = sum(xpos.^2) + sum(zpos.^2);
            dx = bsxfun(@minus,xpos,xpos');
            dz = bsxfun(@minus,zpos,zpos');
            dists = dx.^2 + dz.^2;

            hamiltonian(plot_count) = -omega/2*ang_moment(plot_count) + gval/(8.*pi)*sum(sum(log(dists + eye(Nvorts))));
            radial_densities(plot_count, :) = density_check(xpos,zpos,r0,lscl,ndists);
            %test_dif(plot_count) = norm(u-uref);
            %xtrack = [xtrack;xpos];
            %ztrack = [ztrack;zpos];
            plot_count = plot_count + 1;
            
                
        end                 
            
    end
    
    toc
    
    figure(1)
    %plot(times,ang_moment,'k-',times,ang_moment_ref,'r--','LineWidth',2)
    plot(times,ang_moment,'r--','LineWidth',2)
    

    figure(2)
    %plot(times,ang_moment,'k-',times,ang_moment_ref,'r--','LineWidth',2)
    plot(times,hamiltonian,'r--','LineWidth',2)
    
    %figure(2)
    %plot(times,test_dif,'k-','LineWidth',2)
    
    figure(3)
    Bendixson(xpos,zpos,gval,r0,n_bdry,Nvorts,markersize)

    figure(4)
    surf(10:ndists,times,radial_densities(:,10:end),'LineStyle','none')


end