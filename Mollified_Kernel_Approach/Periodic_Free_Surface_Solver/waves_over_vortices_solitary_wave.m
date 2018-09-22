function waves_over_vortices_solitary_wave(Nx,K,modu,kap,mu,gam,omega,tf)
    
    if modu ~= 0
        Kofk = ellipke(modu);
        Mx = Kofk/kap;     
    else
        Mx = 4;
    end
    disp(Mx)
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    a = 1;
    b = 1;
    sig = 0;
    zoffc = .35;
    av = .5*gam*min(1-zoffc,zoffc);
    %bv = b*mu*gam;
    F = pi*omega*av^2/(4*gam);
    [xpos,zpos,gvals,ep,Nvorts] = initializer(Nx,gam,av,omega);
    zpos = zoffc + zpos/gam;
    
    simul_plot = 0; % Plot during computation.       0 - off, 1 - on
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

        
    uvals = kap*(Xmesh+.5*Mx);
    [~,cn,~] = ellipj(uvals,modu);
    %ceff = 1 + 2/3*mu*kap^2*(2*modu^2-1);
    elipmesh = linspace(0,2*kap,KT+1);
    [~,cnm,~] = ellipj(elipmesh,modu);
    %q0 = -modu^2*kap*(2*kap/KT)*sum((cnm(1:KT)).^2);    
    q0 = 0;
    Q = fft(q0 + 8*modu^2*kap^2*cn.^2);
    eta = Q;    
    
    %{    
    eta = zeros(KT,1);
    phix = init_cond(Xmesh,gam,xpos,zpos,gvals);
    Q = -F*fft(phix);
    Q(Kc:Kuc) = 0;
    %}
    no_dno_term = 40;
    G0 = real(ifft(L1.*Q));
    q = real(ifft([0;-1i*Mx/pi*(1./[1:K -K+1:-1])'.*Q(2:KT)]));
    
    dnonl = dno_maker(real(ifft(eta)),real(ifft(Q)),G0,L1,gam,mu,Kmesh/Mx,no_dno_term);
    k_energy_base = 2*Mx/KT*sum( q.*(G0+dnonl) );
    p_energy_base = 2*Mx/KT*sum( real(ifft(eta)).^2 );
    energy_base = k_energy_base + p_energy_base;
                   
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
    
    samp = 6;
    scnt = 1;
    comx = zeros(floor(nmax/samp),1);
    comt = zeros(floor(nmax/samp),1);
    
    u = [eta;Q;xpos;zpos]; %velocity vector field
    
    % Make folder
    S = make_folder(Nx/2,Nx,K,mu,gam,F,tf,modu,kap);
    clf
    
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %

    if simul_plot
        figure(1)
        Bendixson(Xmesh,Mx,eta,xpos,zpos,gvals,Nvorts,markersize);
        
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
              
        [u,xpos,zpos] = vort_update_on_molly_fourier(Xmesh,gam,mu,ep,u,gvals,L1,no_dno_term,Nvorts,Ehdt,Edt,xpos,zpos,dt,sig,Mx,2*KT);
       
        if(mod(jj,inter)==0)
            eta = u(1:KT);
            Q = u(KT+1:2*KT);
            
            eta = real(ifft(eta));
            G0 = real(ifft(L1.*Q));
            q = real(ifft([0;-1i*Mx/pi*(1./[1:K -K+1:-1])'.*Q(2:KT)]));
            Q = real(ifft(Q));
            
            dnonl = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh/Mx,no_dno_term);
            k_energy_plot(plot_count) = 2*Mx/KT*sum( q.*(G0+dnonl) );
            p_energy_plot(plot_count) = 2*Mx/KT*sum( eta.^2 );
            if modu ~=0
                energy_plot(plot_count) = (k_energy_plot(plot_count) + p_energy_plot(plot_count) - energy_base)/energy_base;
            else
                energy_plot(plot_count) = (k_energy_plot(plot_count) + p_energy_plot(plot_count));
            end
            tm_track(plot_count) = eta(1);
            times(plot_count) = (jj-1)*dt;
            
            plot_count = plot_count + 1;
            
            eta_plot(:,plot_count) = mu*eta;
            
            xtrack = [xtrack;xpos];
            ztrack = [ztrack;zpos];
            gtrack = [gtrack;gvals];
            Vcnt(plot_count) = Nvorts;            
            
            if simul_plot
                Bendixson(Xmesh,Mx,eta,xpos,zpos,gvals,markersize);
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
        
        if(mod(jj,samp)==0)        
            [xpos,zpos,gvals] = recircer_bndry(gvals,xpos,zpos,Nx,Mx);
            Nvorts = length(xpos);
            disp('Current number of vortices is')
            disp(Nvorts)
            u(2*KT+1:2*KT+Nvorts) = xpos;
            u(2*KT+Nvorts+1:2*KT+2*Nvorts) = zpos;
            
            comx(scnt) = com_comp(ep,gam,xpos,zpos,gvals);
            comt(scnt) = jj*dt;
            scnt = scnt + 1;
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
    disp('DNO Resolution')
    disp(norm(dnofin-dnofinn1)/norm(Q))
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Animate waves over vortices
    if modu ~= 0 && ~simul_plot
        etanv = wave_maker_kdv(K,modu,kap,mu,gam,tf);
        pspecnv = log10(abs(fftshift(fft(etanv)))/KT);    

        figure(1)
        gif_my_gif(Xmesh,Mx,eta_plot,mu*etanv,xtrack,ztrack,gtrack,Vcnt,plot_count,S,markersize);
        
        figure(2)
        plot(Xmesh,eta,'k-',Xmesh,etanv,'k--','LineWidth',2)
        h = set(gca,'FontSize',30);
        set(h,'Interpreter','LaTeX')
        xlabel('$x$','Interpreter','LaTeX','FontSize',30)
        ylabel('$\eta(x,t_{f})$','Interpreter','LaTeX','FontSize',30)
        savefig(strcat(S, '/', 'profiles'))    
    end 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % Animate waves over vortices
    figure(2)
    if modu ~= 0
       plot(Xmesh,eta,'k-',Xmesh,etanv,'k--','LineWidth',2)
    else
       plot(Xmesh,eta,'k-','LineWidth',2)
    end
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$x$','Interpreter','LaTeX','FontSize',30)
    ylabel('$\eta(x,t_{f})$','Interpreter','LaTeX','FontSize',30)
    
    %{
    figure(3)
    plot(times,mu*(tm_track-mean(tm_track)),'k','LineWidth',2)
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$t$','Interpreter','LaTeX','FontSize',30)
    ylabel('$\eta(0,t)$','Interpreter','LaTeX','FontSize',30)
    %}
    
    % Plot the power spectrum
    figure(3)
    %plot_pspec(K,S,fftshift(pspec),pspecnv);
    plot(-K+1:K,fftshift(pspec),'k-','LineWidth',2)
    
    % Plot the surface energy
    figure(4)
    plot_energy(plot_count,S,times,energy_plot);    
    
    %{
    tspec = 4*pi*dt*mu^2*(abs(fft(tm_track-mean(tm_track)))/sqrt(no_of_evals)).^2;
    if mod(no_of_evals,2)==0
        tsaxis = 0:no_of_evals/2-1;
    else
        tsaxis = 0:(no_of_evals-1)/2;
    end
    figure(6)
    plot(tsaxis/tf,tspec(1:length(tsaxis)),'k','LineWidth',2)    
    h = set(gca,'FontSize',30);
    set(h,'Interpreter','LaTeX')
    xlabel('$f$','Interpreter','LaTeX','FontSize',30)
    ylabel('$E(f)$','Interpreter','LaTeX','FontSize',30)
    %}
    
    figure(5)
    plot_vorticity(S,ep,gam,xpos,zpos,gvals)
    
    
    figure(6)
    plot_com(S,comt,comx)
    
    disp('Mean relative energy transfer')
    disp(1/(2*(length(energy_plot(1:end-1))))*(energy_plot(1)+energy_plot(end-1)+2*sum(energy_plot(2:end-2))))
    
end