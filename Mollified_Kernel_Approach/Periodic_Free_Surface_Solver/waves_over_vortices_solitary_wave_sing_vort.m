function waves_over_vortices_solitary_wave_sing_vort(K,modu,kap,mu,gam,omega,tf)
    
    if modu ~= 0
        Kofk = ellipke(modu);
        Mx = Kofk/kap;     
    else
        Mx = 8;
    end
    disp(Mx)
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    sig = 0;
    zoffc = .35;
    av = .5*gam*min(1-zoffc,zoffc);
    F = pi*omega*av^2/(4*gam);
    xpos = 0;
    zpos = 0;
    gval = F;
    disp('Froude Number is')
    disp(F)
    zpos = zoffc;
    
    simul_plot = 0; % Plot during computation.       0 - off, 1 - on
    
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
    elipmesh = linspace(0,2*kap,KT+1);
    [~,cnm,~] = ellipj(elipmesh,modu);
    q0 = 0;
    Q = fft(q0 + 8*modu^2*kap^2*cn.^2);
    eta = Q;    
    eta(Kc:Kuc) = 0;
    Q(Kc:Kuc) = 0;
    
    no_dno_term = 20;
    G0 = real(ifft(L1.*Q));
    q = real(ifft([0;-1i*Mx/pi*(1./[1:K -K+1:-1])'.*Q(2:KT)]));
    
    dnonl = dno_maker(real(ifft(eta)),real(ifft(Q)),G0,L1,gam,mu,Kmesh/Mx,no_dno_term);
    k_energy_base = 2*Mx/KT*sum( q.*(G0+dnonl) );
    p_energy_base = 2*Mx/KT*sum( real(ifft(eta)).^2 );
    energy_base = k_energy_base + p_energy_base;
                   
    xtrack = xpos;
    ztrack = zpos;
    
    inter = 1;
    plot_count = 1;
    no_of_evals = round(nmax/inter);
    eta_plot = zeros(KT,no_of_evals+1);
    energy_plot = zeros(no_of_evals,1);
    k_energy_plot = zeros(no_of_evals,1);
    p_energy_plot = zeros(no_of_evals,1);
    times = zeros(no_of_evals,1);
    tm_track = zeros(no_of_evals,1);
    
    comt = zeros(nmax+1,1);
    
    u = [eta;Q;xpos;zpos]; %velocity vector field
    
    % Make folder
    S = make_folder(1,1,K,mu,gam,F,tf,modu,kap,zoffc);
    clf
    
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
    
    tic
    for jj=1:nmax
        % Break if any vortex has a z-value outside of (0,1)
        if(max(zpos) >= 1 || min(zpos) <= 0)
            disp('Out of Bounds!')
            break;
        end
        
        % Now update the vortex positions                                   KTT must be 2*KT because of periodicity!!
              
        [u,xpos,zpos] = vort_update_sing(Xmesh,gam,mu,u,gval,L1,no_dno_term,Ehdt,Edt,xpos,zpos,dt,sig,Mx,2*KT);
       
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
            comt(jj+1) = jj*dt;
            % Stopping Criterion
            pspec = log10(abs(u(1:KT))/KT);            
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
    if modu ~= 0 
        etanv = wave_maker_kdv(K,modu,kap,mu,gam,tf);
        pspecnv = log10(abs(fftshift(fft(etanv)))/KT);    

        figure(1)
        plot(Xmesh,eta,'k-',Xmesh,etanv,'k--','LineWidth',2)
        h = set(gca,'FontSize',30);
        set(h,'Interpreter','LaTeX')
        xlabel('$x$','Interpreter','LaTeX','FontSize',30)
        ylabel('$\eta(x,t_{f})$','Interpreter','LaTeX','FontSize',30)
        %savefig(strcat(S, '/', 'profiles'))    
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
    
    % Plot the power spectrum
    figure(3)
    plot(-K+1:K,fftshift(pspec),'k-','LineWidth',2)
    
    % Plot the surface energy
    figure(4)
    plot_energy(plot_count,S,times,energy_plot);    
    
    figure(5)
    plot_com(S,comt,xtrack)
    
    disp('Mean relative energy transfer')
    disp(1/(2*(length(energy_plot(1:end-1))))*(energy_plot(1)+energy_plot(end-1)+2*sum(energy_plot(2:end-2))))
    
end