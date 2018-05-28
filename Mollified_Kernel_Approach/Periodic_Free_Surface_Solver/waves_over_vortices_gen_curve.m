function waves_over_vortices_gen_curve(Mx,Nx,K,mu,gam,omega,tf)
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    a = 2;
    b = 1;
    av = a*gam*mu;
    bv = b*gam*mu;
    zoffc = .35;
    F = pi*omega*av*bv/gam;
    [xpos,zpos,gvals,ep,Nvorts] = initializer(Nx,gam,av,bv,omega,zoffc);
    
    simul_plot = 0; % Plot during computation.       0 - off, 1 - on
    n_bdry = 0;     % Number of points in cicular boundary.
    stop_crit = 1;  % Stopping criterion.            0 - off, 1 - on
    markersize = 10;
    
    KT = 2*K; %remember periodicity
    % Find the wave numbers to implement the 2/3 de-aliasing throughout
    Kc = floor(2*K/3);
    Kuc = 2*K-Kc+1;
    Kc = Kc+1;
    
    Xmesh = linspace(-Mx,Mx,KT+1);
    Xmesh = Xmesh(1:KT)';
    
    Kmesh = [0:K-1 0 -K+1:-1]';
       
    dt = .05;
    nmax = round(tf/dt);
    
    L1 = -1i*tanh(pi.*gam.*Kmesh/Mx)./gam;
    L2 = -1i*pi*Kmesh/Mx;
    
    [Edt,Ehdt] = get_stuff(Kmesh/Mx,gam,K,KT,dt);
    
    % Here we put in a flat surface and zero background velocity potential
    % Build the quiescent initial velocity potential
    %{
    eta = zeros(KT,1);
    phix = init_cond(Xmesh,gam,xpos,zpos,gvals);
    Q = -F*fft(phix);
    Q(Kc:Kuc) = 0;
    %}    
    % Here we put in the cnoidal solution for KdV in the appropriate
    % coordinates.  
    
    modu = .999;    
    Mval = 1;
    kap = ellipke(modu)*Mval;
    uvals = kap*Xmesh;
    [~,cn,~] = ellipj(uvals,modu);
    ceff = 1 + 2/3*mu*kap^2*(2*modu^2-1);
    %q0 = -modu^2*kap*(2*kap/KT)*sum((cnm(1:KT)).^2);    
    q0 = 0;
    Q = fft(2*q0/3 + 4/3*modu^2*kap^2*cn.^2);
    eta = ceff*Q;    
    plot(Xmesh,real(ifft(eta)),'k','LineWidth',2)
    pause
    eta0 = log10(fftshift(abs(eta))/KT);
    %Q0 = Q;
        
    %{
    amag = 1;
    sig = 2;
    mask = ones(KT,1);
    mask(Kc:Kuc) = 0;
    k0 = 2;
    phi1 = randn(K-1,1);
    phi2 = randn(K-1,1);
    envlpe = amag*exp(-(Kmesh-k0).^2/(2*sig^2)).*mask;
    eta = envlpe.*exp(1i.*[0;phi1;0;conj(flipud(phi1))]);
    Q = envlpe.*exp(1i.*[0;phi2;0;conj(flipud(phi2))]);
    eta0 = log10(fftshift(eta)/KT);
    %}    
    xtrack = xpos;
    ztrack = zpos;
    gtrack = gvals;
    
    inter = 1;
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
    
    no_dno_term = 30;
    
    u = [eta;Q;xpos;zpos]; %velocity vector field
    
    % Make folder
    S = make_folder(Nx/2,Nx,K,mu,gam,F,tf,a/b);
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
              
        [u,xpos,zpos] = vort_update_on_molly_fourier(Xmesh,gam,mu,ep,u,gvals,L1,no_dno_term,Nvorts,Ehdt,Edt,xpos,zpos,dt,Mx,2*KT);
       
        if(mod(jj,inter)==0)
            eta = u(1:KT);
            Q = u(KT+1:2*KT);
            
            eta = real(ifft(eta));
            G0 = real(ifft(L1.*Q));
            q = real(ifft([0;-1i*Mx/pi*(1./[1:K -K+1:-1])'.*Q(2:KT)]));
            Q = real(ifft(Q));
            
            dnonl = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh/Mx,no_dno_term);
            k_energy_plot(plot_count) = Mx/KT*sum( q.*(G0+dnonl) );
            p_energy_plot(plot_count) = Mx/KT*sum( eta.^2 );
            energy_plot(plot_count) = k_energy_plot(plot_count) + p_energy_plot(plot_count);
            tm_track(plot_count) = eta(1);
            times(plot_count) = (jj-1)*dt;
            
            plot_count = plot_count + 1;
            
            eta_plot(:,plot_count) = mu*eta;
            
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
            %{
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
            %}
        end
        
        if(mod(jj,8*inter)==0)        
            [xpos,zpos,gvals] = recircer_bndry(gvals,xpos,zpos,Nx);
            Nvorts = length(xpos);
            disp('Current number of vortices is')
            disp(Nvorts)
            u(2*KT+1:2*KT+Nvorts) = xpos;
            u(2*KT+Nvorts+1:2*KT+2*Nvorts) = zpos;
        end
        
    end
    toc
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    Q = u(KT+1:2*KT);
    
    eta = real(ifft(u(1:KT)));
    G0 = real(ifft(L1.*Q));
    Q = real(ifft(Q));
    
    dnofin = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh/Mx,no_dno_term);
    dnofinn1 = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh/Mx,no_dno_term-1);
    disp(norm(dnofin-dnofinn1)/norm(Q))
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Animate waves over vortices
    if ~simul_plot
        figure(1)
        gif_my_gif(Xmesh,eta_plot,xtrack,ztrack,gtrack,n_bdry,Vcnt,plot_count,S,markersize);
    end

    figure(2)
    plot(times,mu*(tm_track-mean(tm_track)),'k','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$t$','Interpreter','LaTeX','FontSize',30)
    ylabel('$\eta(0,t)$','Interpreter','LaTeX','FontSize',30)

    % Plot the power spectrum
    figure(3)
    plot_pspec(K,S,fftshift(pspec),eta0);
    
    % Plot the surface energy
    figure(4)
    plot_energy(plot_count,S,times,energy_plot,p_energy_plot,k_energy_plot);
    
    tspec = 4*pi*dt*mu^2*(abs(fft(tm_track-mean(tm_track)))/sqrt(no_of_evals)).^2;
    if mod(no_of_evals,2)==0
        tsaxis = 0:no_of_evals/2-1;
    else
        tsaxis = 0:(no_of_evals-1)/2;
    end
    figure(5)
    plot(tsaxis/tf,tspec(1:length(tsaxis)),'k','LineWidth',2)    
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$f$','Interpreter','LaTeX','FontSize',30)
    ylabel('$E(f)$','Interpreter','LaTeX','FontSize',30)

end