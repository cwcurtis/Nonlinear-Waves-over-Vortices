function waves_over_vortices_on_molly_fourier(rows,cols,K,mu,gam,ep,F,tf,Ntrunc)
    close all
    Gamma = 4;
    xpos1 = 0; % initial x position of vortex formation
    zpos1 = .25; % initial z position of vortex formation
    z_range = 4; % vertical range or diameter of vortex formation
    x_range = 4; % horizontal range of vortex formation (if sqr_grid == 0)
    
    coords = 'polar';
%    coords = 'cartesian';
%    p_shape = 'rose';
%    p_shape = 'cardioid';
    p_shape = 'circle';
%    pattern = 'gradient';
%    pattern = 'clustered';
%    pattern = 'checkerboard';
%    pattern = 'ver_stripes';
%    pattern = 'hor_stripes';
    pattern = 'uniform';
    simul_plot = 0; % Plot during computation.       0 - off, 1 - on
    sqr_grid = 1;   % Same grid spacing for x and z. 0 - off, 1 - on
    disk = 0;       % Disk shape for cartesian grid. 0 - off, 1 - on
    g2 = 0;         % Secondary grid.                0 - off, 1 - on
    auto_gvals = 1; % Automatically assign gvals.    0 - off, 1 - on
    stop_crit = 1;  % Stopping criterion.            0 - off, 1 - on
    n_bdry = 0;     % Number of points in cicular boundary.
    markersize = 10;
    tic
    
    % rows - number of rows in vortex set
    % cols - number of columns in vortex set
    % K - number of spectral modes 
    % mu - the parameter mu from the write-up
    % gam - the parameter gamma from the write-up
    % tf - final time to which to run the simulation
    
    % This program calls several sub-functions.  These are:
    
    % dno_maker - builds the Dirichlet-to-Neumann (DNO) expansion 
    
    % force_terms - builds Pv, Ev, and the nonlinearity associated with
    % Bernoulli's equation
    
    % vort_update - used to run the Runge-Kutta scheme for updating the
    % vortex positions
    
    % Movie_Maker_1d - used to produce animations of the surface profile.

    KT = 2*K; %remember periodicity
    % Find the wave numbers to implement the 2/3 de-aliasing throughout
    Kc = floor(2*K/3);
    Kuc = 2*K-Kc+1;
    Kc = Kc+1;
    
    Xmesh = linspace(-1,1,KT+1);
    Xmesh = Xmesh(1:KT)';
    
    Kmesh = [0:K-1 0 -K+1:-1]';
    
    r_in = 0.8; % Ratio of the radii of the two grids
    d_in = 0; % Ratio of the densities of the two grids
    rows_in = ceil(d_in*rows);
    cols_in = ceil(d_in*cols);
    Nvorts = rows*cols+rows_in*cols_in+n_bdry;
    c1 = Gamma/(rows*cols);
    c2 = -c1;
    c3 = c1;
    c4 = -c1;
    c5 = 4*c1;
    
    if auto_gvals == 1
        % Automatically set vortex strengths
        gvals = set_gvals(rows,cols,rows_in,cols_in,Nvorts,c1,c2,c3,c4,c5,pattern,g2);
    else
        % Manually set vortex strengths
        gvals(1) = c1;
        gvals(3) = -c1;
        gvals(5) = -c1;
        gvals(7) = c1;
        gvals(9) = -c1;
        gvals(13) = c1;
    end
    
    % Choose grid spacings
    if strcmp(coords,'cartesian')
        dz = z_range*mu*gam/max(rows-1,1);
    else
        dz = z_range*mu*gam/rows;
    end
    
    if rows == 1 || sqr_grid == 0
        dx = x_range*mu*gam/max(cols-1,1);
    else
        dx = dz;
    end

    % Initialize positions
    [xpos,zpos] = init_pos(rows,cols,rows_in,cols_in,n_bdry,Nvorts,xpos1,zpos1,dx,dz,coords,g2,p_shape);
    
    % Disk Shape
    if disk == 1
        for ii = 1:Nvorts
            r_i = (2*(xpos(ii)-xpos1)/(cols*dx))^2+(2*(zpos(ii)-zpos1)/(rows*dz))^2;
            gvals(ii) = gvals(ii)*(r_i < 1);
            
            % Secondary Disk
            if g2 == 1
                gvals(ii) = gvals(ii)*(((r_i >= r_in) && (ii <= rows*cols)) || ((r_i < r_in) && (ii > rows*cols)));
            end
        end
    end
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    
    dt = 1e-2;
    nmax = round(tf/dt);
    
    L1 = -1i*tanh(pi.*gam.*Kmesh)./gam;
    L2 = -1i*pi*Kmesh;
    
    [Edt,Ehdt] = get_stuff(Kmesh,gam,K,KT,dt);
    
    % Here we put in a flat surface and zero background velocity potential
    
    % Build the quiescent initial velocity potential
    eta = zeros(KT,1);
    phix = init_cond(Xmesh,gam,xpos,zpos,gvals);
    Q = -F*fft(phix);
    Q(Kc:Kuc) = 0;
    Q0 = Q;
    
    % Here we put in the cnoidal solution for KdV in the appropriate
    % coordinates.
    %{
    modu = .9;
    Mval = 2;
    kap = ellipke(modu)/Mval;
    uvals = kap*Xmesh;
    [~,cn,~] = ellipj(uvals,modu);
    
    ceff = 1 + 2/3*mu*kap^2*(2*modu^2-1);
    
    elipmesh = linspace(0,2*kap,KT+1);
    [~,cnm,~] = ellipj(elipmesh,modu);
    q0 = -modu^2*kap*(2*kap/KT)*sum((cnm(1:KT)).^2);
    
    %q0 = 0;
    Q = fft(2*q0/3 + 4/3*modu^2*kap^2*cn.^2);
    
    eta = ceff*Q;
    %}
    eta0 = log10(fftshift(abs(eta))/KT);
    eta0p = real(ifft(eta));
    
    xtrack = zeros(Nvorts,nmax+1);
    ztrack = zeros(Nvorts,nmax+1);
    
    xtrack(:,1) = xpos;
    ztrack(:,1) = zpos;
    
    inter = 10;
    plot_count = 1;
    no_of_evals = round(nmax/inter);
    eta_plot = zeros(KT,no_of_evals+1);
    energy_plot = zeros(no_of_evals,1);
    k_energy_plot = zeros(no_of_evals,1);
    p_energy_plot = zeros(no_of_evals,1);
    times = zeros(no_of_evals,1);
    
    no_dno_term = 15;
    
    u = [eta;Q;xpos;zpos]; %velocity vector field
    
    % Make folder
    S = make_folder(rows,cols,K,mu,gam,F,tf,Gamma,coords,pattern);
    filename = strcat(S, '/', '/waves_over_vortices.gif');
    clf
    
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %

    if simul_plot
        figure(1)
        Bendixson(Xmesh,eta,xpos,zpos,gvals,n_bdry,Nvorts,markersize);
        
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','DelayTime',0,'loopcount',inf);
    end
    
    tic
    for jj=1:nmax
        % Break if any vortex has a z-value outside of (0,1)
        if(max(zpos) >= 1 || min(zpos) <= 0)
            break;
        end
        
        % Now update the vortex positions                                   KTT must be 2*KT because of periodicity!!
        [u,xpos,zpos] = vort_update_on_molly_fourier(Xmesh,gam,mu,ep,F,u,gvals,L1,no_dno_term,Nvorts,Ntrunc,Ehdt,Edt,xpos,zpos,dt,2*KT);
        
        xtrack(:,jj+1) = xpos;
        ztrack(:,jj+1) = zpos;
        
        if(mod(jj,inter)==0)
            eta = u(1:KT);
            Q = u(KT+1:2*KT);
            
            eta = real(ifft(eta));
            G0 = real(ifft(L1.*Q));
            q = real(ifft([0;-1i/pi*(1./[1:K -K+1:-1])'.*Q(2:KT)]));
            Q = real(ifft(Q));
            
            dnonl = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh,no_dno_term);
            k_energy_plot(plot_count) = 1/KT*sum( q.*(G0+dnonl) );
            p_energy_plot(plot_count) = 1/KT*sum( eta.^2 );
            energy_plot(plot_count) = k_energy_plot(plot_count) + p_energy_plot(plot_count);
            times(plot_count) = (jj-1)*dt;
            
            plot_count = plot_count + 1;
            
            eta_plot(:,plot_count) = eta;
            
            if simul_plot
                Bendixson(Xmesh,eta,xpos,zpos,gvals,n_bdry,Nvorts,markersize);

                frame = getframe(1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                imwrite(imind,cm,filename,'gif','DelayTime',0,'writemode','append');
            end
            
            % Stopping Criterion
            pspec = log10(fftshift(abs(u(1:KT)))/KT);
            if stop_crit == 1
                pspecfit = polyfit((1:floor(K*2/3-1))',pspec(K+2:floor(K*5/3)),1);
                if (pspecfit(1) > -.1)
                    break;
                end
            end
        end
    end
    toc
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    Q = u(KT+1:2*KT);
    
    eta = real(ifft(u(1:KT)));
    G0 = real(ifft(L1.*Q));
    Q = real(ifft(Q));
    
    dnofin = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh,no_dno_term);
    dnofinn1 = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh,no_dno_term-1);
    disp(norm(dnofin-dnofinn1)/norm(Q))
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Animate waves over vortices
    if ~simul_plot
        figure(1)
        gif_my_gif(Xmesh,eta_plot,xtrack,ztrack,gvals,n_bdry,Nvorts,inter,plot_count,filename,markersize);
    end

    % Plot the power spectrum
    figure(2)
    plot_pspec(K,S,pspec,eta0);
    
    % Plot the surface energy
    figure(3)
    plot_energy(plot_count,S,times,energy_plot,p_energy_plot,k_energy_plot);
    
    toc
end