function waves_over_vortices_gen_curve(Nx,K,mu,gam,omega,ep,tf,Ntrunc)
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    cfun = @(x,z,gv,av,bv,zoff) (x.^2/av^2 + gv^2*(z-zoff).^2/bv^2);
    avc = 1/25;
    bvc = 1/25;
    zoffc = .35;
    F = pi*omega*avc*bvc/gam;
    [xpos,zpos,gvals,Gamma,Nvorts] = initializer(Nx,F,cfun,gam,avc,bvc,zoffc);
    
    simul_plot = 0; % Plot during computation.       0 - off, 1 - on
    n_bdry = 0;     % Number of points in cicular boundary.
    stop_crit = 1;  % Stopping criterion.            0 - off, 1 - on
    markersize = 10;
    
    KT = 2*K; %remember periodicity
    % Find the wave numbers to implement the 2/3 de-aliasing throughout
    Kc = floor(2*K/3);
    Kuc = 2*K-Kc+1;
    Kc = Kc+1;
    
    Xmesh = linspace(-1,1,KT+1);
    Xmesh = Xmesh(1:KT)';
    
    Kmesh = [0:K-1 0 -K+1:-1]';
       
    dt = .02;
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
    
    eta0 = log10(fftshift(abs(eta))/KT);
    eta0p = real(ifft(eta));
    
    xtrack = xpos;
    ztrack = zpos;
    gtrack = gvals;
    
    inter = 10;
    plot_count = 1;
    no_of_evals = round(nmax/inter);
    eta_plot = zeros(KT,no_of_evals+1);
    energy_plot = zeros(no_of_evals,1);
    k_energy_plot = zeros(no_of_evals,1);
    p_energy_plot = zeros(no_of_evals,1);
    times = zeros(no_of_evals,1);
    tm_track = zeros(no_of_evals,1);
    Vcnt = zeros(no_of_evals+1,1);
    Vcnt(1) = Nvorts;
    
    no_dno_term = 15;
    
    u = [eta;Q;xpos;zpos]; %velocity vector field
    
    % Make folder
    S = make_folder(Nx/2,Nx,K,mu,gam,F,tf,Gamma);
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
            disp('Out of Bounds!')
            break;
        end
        
        % Now update the vortex positions                                   KTT must be 2*KT because of periodicity!!
              
        [u,xpos,zpos] = vort_update_on_molly_fourier(Xmesh,gam,mu,ep,F,u,gvals,L1,no_dno_term,Nvorts,Ntrunc,Ehdt,Edt,xpos,zpos,dt,2*KT);
       
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
            tm_track(plot_count) = eta(1);
            times(plot_count) = (jj-1)*dt;
            
            plot_count = plot_count + 1;
            
            eta_plot(:,plot_count) = eta;
            
            xtrack = [xtrack;xpos];
            ztrack = [ztrack;zpos];
            gtrack = [gtrack;gvals];
            Vcnt(plot_count) = Nvorts;            
            
            if simul_plot
                Bendixson(Xmesh,eta,xpos,zpos,gvals,n_bdry,Nvorts,markersize);
                frame = getframe(1);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                imwrite(imind,cm,filename,'gif','DelayTime',0,'writemode','append');
            end
            
            % Stopping Criterion
            pspec = log10(abs(u(1:KT))/KT);
            if stop_crit == 1
                kvals = (11:floor(Kc))';
                ntot = length(kvals);
                fvals = pspec(11:floor(Kc));
                slp = (ntot*sum(kvals.*fvals)-sum(kvals)*sum(fvals))/(ntot*sum(kvals.^2)-(sum(kvals))^2);
                if (slp > -.1)
                    disp('Too Bad Your Spectrum is Mad')
                    break;
                end
            end            
             
        end
        %{
        if(mod(jj,inter)==0)        
            [xpos,zpos,gvals] = recircer(gvals,xpos,zpos,Nx);
            Nvorts = length(xpos);
            disp('Current number of vortices is')
            disp(Nvorts)
            u(2*KT+1:2*KT+Nvorts) = xpos;
            u(2*KT+Nvorts+1:2*KT+2*Nvorts) = zpos;
        end
        %}
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
        gif_my_gif(Xmesh,eta_plot,xtrack,ztrack,gtrack,n_bdry,Vcnt,plot_count,filename,markersize);
    end

    % Plot the power spectrum
    figure(2)
    plot_pspec(K,S,fftshift(pspec),eta0);
    
    % Plot the surface energy
    figure(3)
    plot_energy(plot_count,S,times,energy_plot,p_energy_plot,k_energy_plot);
    
    tspec = 4*pi*tf*(abs(fft(tm_track))/sqrt(no_of_evals)).^2;
    if mod(no_of_evals,2)==0
        tsaxis = 0:no_of_evals/2-1;
    else
        tsaxis = 0:(no_of_evals-1)/2;
    end
    figure(4)
    plot(tsaxis/tf,log10(tspec(1:length(tsaxis))),'k','LineWidth',2)    
end