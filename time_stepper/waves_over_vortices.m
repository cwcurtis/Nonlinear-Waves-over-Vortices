function waves_over_vortices(K,mu,gam,F,tf)

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

    KT = 2*K;
    % Find the wave numbers to implement the 2/3 de-aliasing throughout
    Kc = floor(2*K/3);
    Kuc = 2*K-Kc+1;
    Kc = Kc+1;
    
    Xmesh = linspace(-1,1,KT+1);
    Xmesh = Xmesh(1:KT)';    
    
    Kmesh = [0:K-1 0 -K+1:-1]';
    
    % Necessary information for running the two vortex case 
    
    Nvorts = 2;
    c1 = 1;
    c2 = -c1;
    gvals = [c1;c2];
    
    xpos = [-mu*gam;mu*gam];
    
    zpos1 = .25;
    zpos = [zpos1;zpos1];
    
    % Necessary information for running the four vortex case 
    %{
    Nvorts = 4;
    c1 = 1;
    c2 = -c1;
    c3 = c1;
    c4 = -c1;
    gvals = [c1;c2;c3;c4];
    
    xpos = [-2*mu*gam-mu*gam/2;-2*mu*gam+mu*gam/2;2*mu*gam-mu*gam/2;2*mu*gam+mu*gam/2];
    
    zpos1 = .25;
    
    zpos = [zpos1;zpos1;zpos1;zpos1];
    %}
    % Necessary information for running the six vortex case 
    %{
    Nvorts = 6;
    c1 = 1;
    c2 = c1;
    c3 = c1;
    c4 = -c1;
    c5 = -c1;
    c6 = -c1;
    gvals = [c1;c2;c3;c4;c5;c6];
    
    xpos = [-mu*gam;-2*mu*gam/3;-mu*gam/3;mu*gam/3;2*mu*gam/3;mu*gam];
    
    zpos1 = .25;
    
    zpos = [zpos1;zpos1-.05;zpos1;zpos1;zpos1-.05;zpos1];
    %}      
    disp([xpos zpos])
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.  
    
    dt = 5e-3;
    nmax = round(tf/dt);
    
    L1 = -1i*tanh(pi.*gam.*Kmesh)./gam;
    L2 = -1i*pi*Kmesh;
    
    omf = @(k) 1i*sqrt(pi.*k.*tanh(pi.*gam.*k)./gam);
    disp(dt*max(abs(omf(Kmesh))))
    vmatf = @(k) sqrt(pi.*gam.*k./tanh(pi.*gam.*k));
    vmatv = sign(Kmesh).*vmatf(Kmesh);
    vmatv(1) = 1;
    vmatv(K+1) = 1;
    
    Vmat = [eye(KT) eye(KT);-diag(vmatv) diag(vmatv)];
    Vmati = 1/2*[eye(KT) -diag(vmatv.^(-1)); eye(KT) diag(vmatv.^(-1))];
    
    Edt = Vmat*diag([exp(dt*omf(Kmesh));exp(-dt*omf(Kmesh))])*Vmati;
    Ehdt = Vmat*diag([exp(dt/2*omf(Kmesh));exp(-dt/2*omf(Kmesh))])*Vmati;
    
    % Here we put in a flat surface and zero background velocity potential
    
    eta = zeros(KT,1);
    
    % Build the quiescent initial velocity potential 
    phix = init_cond(Xmesh,gam,xpos,zpos,gvals);
    Q = -F*fft(phix);
    Q(Kc:Kuc) = 0;
    Q0 = Q;
    
    % Here we put in the cnoidal solution for KdV in the appropriate
    % coordinates.  
    %{
    modu = .3;    
    kap = ellipke(modu);
    uvals = kap*Xmesh;
    [~,cn,~] = ellipj(uvals,modu);
    elipmesh = linspace(0,2*kap,KT+1);
    [~,cnm,~] = ellipj(elipmesh,modu);
    
    ceff = 1 + 2/3*mu*kap^2*(2*modu^2-1);
    
    %q0 = -modu^2*kap*(2*kap/KT)*sum((cnm(1:KT)).^2);
    q0 = 0;
    Q = fft(2*q0/3 + 4/3*modu^2*kap^2*cn.^2);
    
    eta = ceff*Q;
    %}
    eta0 = log10(fftshift(abs(eta))/KT);
    
    xtrack = zeros(Nvorts,nmax+1);
    ztrack = zeros(Nvorts,nmax+1);
    
    xtrack(:,1) = xpos;
    ztrack(:,1) = zpos;
    
    nthrd = round(nmax/3);
    count = 1;
    
    inter = 10;
    plot_count = 1;
    no_of_evals = round(nmax/inter);
    energy_plot = zeros(no_of_evals,1); 
    times = zeros(no_of_evals,1);
    no_dno_term = 17;
    
    u = [eta;Q;xpos;zpos];
    
    for jj=1:nmax
            
        % Now update the vortex positions
        k1 = dt*force_terms(Xmesh,gam,mu,F,u,gvals,L1,no_dno_term,Nvorts);
        uh1 = [Ehdt*(u(1:2*KT)+k1(1:2*KT)/2);xpos+k1(2*KT+1:2*KT+Nvorts)/2;zpos+k1(2*KT+Nvorts+1:2*KT+2*Nvorts)/2];
        
        k2 = dt*force_terms(Xmesh,gam,mu,F,uh1,gvals,L1,no_dno_term,Nvorts);
        uh2 = [(Ehdt*u(1:2*KT)+k2(1:2*KT)/2);xpos+k2(2*KT+1:2*KT+Nvorts)/2;zpos+k2(2*KT+Nvorts+1:2*KT+2*Nvorts)/2];
        
        k3 = dt*force_terms(Xmesh,gam,mu,F,uh2,gvals,L1,no_dno_term,Nvorts);
        uh3 = [(Edt*u(1:2*KT)+Ehdt*k3(1:2*KT));xpos+k3(2*KT+1:2*KT+Nvorts);zpos+k3(2*KT+Nvorts+1:2*KT+2*Nvorts)];
        
        k4 = dt*force_terms(Xmesh,gam,mu,F,uh3,gvals,L1,no_dno_term,Nvorts);
        
        u(1:2*KT) = Edt*u(1:2*KT) + (Edt*k1(1:2*KT) + 2*Ehdt*(k2(1:2*KT) + k3(1:2*KT)) + k4(1:2*KT))/6;
        u(2*KT+1:2*KT+2*Nvorts) = u(2*KT+1:2*KT+2*Nvorts) + ( k1(2*KT+1:2*KT+2*Nvorts) + 2*(k2(2*KT+1:2*KT+2*Nvorts)+k3(2*KT+1:2*KT+2*Nvorts)) + k4(2*KT+1:2*KT+2*Nvorts) )/6;
        
        xpos = u(2*KT+1:2*KT+Nvorts);
        zpos = u(2*KT+Nvorts+1:2*KT+2*Nvorts);
        
        if(max(zpos) >= 1)
            break;
        else
                
            indsp = find(xpos>1);
            xpos(indsp) = xpos(indsp) - 1;
        
            indsn = find(xpos<-1);
            xpos(indsn) = 1 + xpos(indsn);
        
            xtrack(:,jj+1) = xpos;
            ztrack(:,jj+1) = zpos;     
        
            % Now update all previous time step values for next iteration
        
            count = count + 1;
           
        end
        
        % Add to the movie plot matrix for the purpose of making animations
        
        if jj==nthrd
            fst_plot = real(ifft(u(1:KT)));            
        end
        
        if jj==2*nthrd
            scd_plot = real(ifft(u(1:KT)));
        end
        
        if(mod(jj,inter)==0)
            eta = u(1:KT);
            Q = u(KT+1:2*KT);
            eta = real(ifft(eta));
            G0 = real(ifft(L1.*Q));        
            q = real(ifft([0;-1i/pi*(1./[1:K -K+1:-1])'.*Q(2:KT)]));
            Q = real(ifft(Q));
            dnonl = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh,no_dno_term);
            energy_plot(plot_count) = 1/KT*sum( q.*(G0+dnonl) + eta.^2 );
            times(plot_count) = (jj-1)*dt;
            plot_count = plot_count + 1;
       end
        
    end
    
    eta = u(1:KT);
    pspec = log10(fftshift(abs(eta))/KT);
    Q = u(KT+1:2*KT);
    
    eta = real(ifft(eta));
    G0 = real(ifft(L1.*Q));        
    Q = real(ifft(Q));
    
    dnofin = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh,no_dno_term);
    dnofinn1 = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh,no_dno_term-1);
    disp(norm(dnofin-dnofinn1)/norm(Q))    
        
    % If you want to compare against the analytic formula for the linear
    % motion
    %{
    linrep = linear_response(K,Q0,gam,F,xpos(2),zpos(2),tf);
    figure(1)
    plot(Xmesh,eta,'r--',Xmesh,linrep,'k.','LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylabel('\eta(x,t)','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    %}
    % Plot the surface profile
    
    clf
    figure(1)
    plot(Xmesh,fst_plot,'k-.',Xmesh,scd_plot,'k--',Xmesh,eta,'k','LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylabel('\eta(x,t)','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    
    % Plot the power spectrum
    
    figure(2)
    plot(-K+1:K,pspec,'k',-K+1:K,eta0,'LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('k','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylabel('\eta(k,t)','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylim([-20 0])
    xlim([-K+1 K])
    
    % Plotting the paths in the two vortex case  
        
    figure(3)
    hold on
    plot(xtrack(1,1:count),ztrack(1,1:count),'k',xtrack(2,1:count),ztrack(2,1:count),'k','LineWidth',2)
    plot(xtrack(1,1),ztrack(1,1),'.','MarkerSize',26','color',[0.8 0.8 0.8]);
    plot(xtrack(2,1),ztrack(2,1),'.','MarkerSize',26','color',[0.8 0.8 0.8]);
    plot(xtrack(1,end),ztrack(1,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
    plot(xtrack(2,end),ztrack(2,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
    hold off
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylabel('z','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    
    % Plotting the paths in the four vortex case   
    %{  
    figure(3)
    hold on
    plot(xtrack(1,1:count),ztrack(1,1:count),'k',xtrack(2,1:count),ztrack(2,1:count),'k-.',...
         xtrack(3,1:count),ztrack(3,1:count),'k-.',xtrack(4,1:count),ztrack(4,1:count),'k','LineWidth',2)
    
    plot(xtrack(1,1),ztrack(1,1),'.','MarkerSize',26','color',[0.8 0.8 0.8]);
    plot(xtrack(2,1),ztrack(2,1),'.','MarkerSize',26','color',[0.8 0.8 0.8]);
    plot(xtrack(3,1),ztrack(3,1),'.','MarkerSize',26','color',[0.8 0.8 0.8]);
    plot(xtrack(4,1),ztrack(4,1),'.','MarkerSize',26','color',[0.8 0.8 0.8]);
    
    plot(xtrack(1,end),ztrack(1,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
    plot(xtrack(2,end),ztrack(2,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
    plot(xtrack(3,end),ztrack(3,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
    plot(xtrack(4,end),ztrack(4,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
    hold off
    
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylabel('z','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    %}    
    % Plotting the paths in the six vortex case
    %{
    figure(3)
    plot(xtrack(1,1:count),ztrack(1,1:count),xtrack(2,1:count),ztrack(2,1:count),xtrack(3,1:count),ztrack(3,1:count),xtrack(4,1:count),ztrack(4,1:count),xtrack(5,1:count),ztrack(5,1:count),xtrack(6,1:count),ztrack(6,1:count),'LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylabel('z','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    %}   
    
    figure(4)
    plot(times,energy_plot,'k','LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('t','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylabel('E(t)','FontName','Helvetica','FontSize',30,'FontWeight','bold')
   
end