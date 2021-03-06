function wave_maker(K,Mx,mu,gam,tf)
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    
    simul_plot = 0; % Plot during computation.       0 - off, 1 - on
    stop_crit = 1;  % Stopping criterion.            0 - off, 1 - on
    
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
    kap = ellipke(modu)*Mval/Mx;
    uvals = kap*Xmesh;
    [~,cn,~] = ellipj(uvals,modu);
    ceff = 1 + 2/3*mu*kap^2*(2*modu^2-1);
    elipmesh = linspace(0,2*kap,KT+1);
    [~,cnm,~] = ellipj(elipmesh,modu);
    %q0 = -modu^2*kap*(2*kap/KT)*sum((cnm(1:KT)).^2);    
    q0 = 0;
    Q = fft(2*q0/3 + 4/3*modu^2*kap^2*cn.^2);
    eta = ceff*Q;    
    eta0 = log10(fftshift(eta)/KT);    
    plot(Xmesh,real(ifft(eta)),'k','LineWidth',2)
    pause
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
    inter = 1;
    plot_count = 1;
    
    no_of_evals = round(nmax/inter);
    eta_plot = zeros(KT,no_of_evals+1);
    energy_plot = zeros(no_of_evals,1);
    k_energy_plot = zeros(no_of_evals,1);
    p_energy_plot = zeros(no_of_evals,1);
    times = zeros(no_of_evals,1);
    tm_track = zeros(no_of_evals,1);
    
    no_dno_term = 15;
    
    u = [eta;Q]; %velocity vector field
    
    % Make folder
    S = make_folder(K,mu,gam,tf);
    clf
    
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %

    if simul_plot
        figure(1)
        Bendixson(Xmesh,Mx,eta);
        
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','DelayTime',0,'loopcount',inf);
    end
    
    tic
    for jj=1:nmax
        
        u = vort_update_on_molly_fourier(Xmesh,Mx,gam,mu,u,L1,no_dno_term,Ehdt,Edt,dt);
       
        if(mod(jj,inter)==0)
            eta = u(1:KT);
            Q = u(KT+1:2*KT);
            
            eta = real(ifft(eta));
            G0 = real(ifft(L1.*Q));
            q = real(ifft([0;-1i*Mx/pi*(1./[1:K -K+1:-1])'.*Q(2:KT)]));
            Q = real(ifft(Q));
            
            dnonl = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh,no_dno_term);
            k_energy_plot(plot_count) = 1/KT*sum( q.*(G0+dnonl) );
            p_energy_plot(plot_count) = 1/KT*sum( eta.^2 );
            energy_plot(plot_count) = k_energy_plot(plot_count) + p_energy_plot(plot_count);
            tm_track(plot_count) = eta(1);
            times(plot_count) = (jj-1)*dt;
            
            plot_count = plot_count + 1;
            
            eta_plot(:,plot_count) = mu*eta;
            
            if simul_plot
                Bendixson(Xmesh,Mx,mu*eta);
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
        gif_my_gif(Xmesh,Mx,eta_plot,plot_count,S);
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
    plot(tsaxis/tf, tspec(1:length(tsaxis)),'k','LineWidth',2)    
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$f$','Interpreter','LaTeX','FontSize',30)
    ylabel('$E(f)$','Interpreter','LaTeX','FontSize',30)

end