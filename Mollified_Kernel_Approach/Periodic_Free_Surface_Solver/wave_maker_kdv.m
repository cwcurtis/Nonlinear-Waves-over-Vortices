function eta = wave_maker_kdv(K,modu,kap,mu,gam,tf)
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    
    %simul_plot = 0; % Plot during computation.       0 - off, 1 - on
    %stop_crit = 1;  % Stopping criterion.            0 - off, 1 - on
    
    KT = 2*K; %remember periodicity
    % Find the wave numbers to implement the 2/3 de-aliasing throughout
    Kc = floor(2*K/3);
    Kuc = 2*K-Kc+1;
    Kc = Kc+1;
    
    Kofk = ellipke(modu);
    Mx = Kofk/kap;     
    disp(Mx)
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
    ceff = 1 + 2/3*mu*kap^2*(2*modu^2-1);
    elipmesh = linspace(0,2*kap,KT+1);
    [~,cnm,~] = ellipj(elipmesh,modu);
    %q0 = -modu^2*kap*(2*kap/KT)*sum((cnm(1:KT)).^2);    
    q0 = 0;
    Q = fft(q0 + 16*modu^2*kap^2*cn.^2)/6^(1/3);
    eta = ceff*Q;    
    
    no_dno_term = 20;
    
    u = [eta;Q]; %velocity vector field    
    
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % %

    for jj=1:nmax
       u = rk4_time_step(Xmesh,Mx,gam,mu,u,L1,no_dno_term,Ehdt,Edt,dt);
    end
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    Q = u(KT+1:2*KT);
    
    eta = real(ifft(u(1:KT)));
    G0 = real(ifft(L1.*Q));
    Q = real(ifft(Q));        
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
end