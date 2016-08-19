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
    Dx = 1i*pi*Kmesh;
    
    % Necessary information for running the four vortex case 
    
    Nvorts = 4;
    c1 = 1;
    c2 = c1;
    c3 = -c1;
    c4 = -c1;
    gvals = [c1;c2;c3;c4];
    
    xpos = [-2*mu*gam-mu*gam/2;-2*mu*gam+mu*gam/2;2*mu*gam-mu*gam/2;2*mu*gam+mu*gam/2];
    
    zpos1 = .25;
    
    zpos = [zpos1;zpos1;zpos1;zpos1];
    
    % Necessary information for running the two vortex case 
    %{
    Nvorts = 2;
    c1 = 1;
    c2 = -c1;
    gvals = [c1;c2];
    
    xpos = [-mu*gam;mu*gam];
    
    zpos1 = .25;
    zpos = [zpos1;zpos1];
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
    modu = .2;    
    kap = ellipke(modu);
    uvals = kap*Xmesh;
    [~,cn,~] = ellipj(uvals,modu);
    
    ceff = 1 + 2/3*mu*kap^2*(2*modu^2-1);
    
    Q = fft(4/3*modu^2*kap^2*cn.^2);
    
    eta = ceff*Q;
    %}
    xtrack = zeros(Nvorts,nmax+1);
    ztrack = zeros(Nvorts,nmax+1);
    
    xtrack(:,1) = xpos;
    ztrack(:,1) = zpos;
    
    inter = 10; % Number of time steps between frames of movie
    no_of_plots = round(nmax/inter); % Number of plots in movie
    movie_plot = zeros(no_of_plots,KT); % Movie storage
    plot_count = 1;
    count = 1;
    
    no_dno_term = 31;
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
        
        if(max(zpos) >= .95)
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
        %{
        if(mod(jj,inter)==0)
            movie_plot(plot_count,:) = real(ifft(u(1:KT)));
            plot_count = plot_count + 1;
        end
        %}
    end
    
    eta = u(1:KT);
    pspec = log10(fftshift(abs(eta))/KT);
    Q = u(KT+1:2*KT);
    
    eta = real(ifft(eta));
    Q = real(ifft(Q));
    G0 = real(ifft(L1.*Q));        
    
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
    
    figure(1)
    plot(Xmesh,eta,'k','LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylabel('\eta(x,t)','FontName','Helvetica','FontSize',30,'FontWeight','bold')
       
    % Plot the power spectrum
    
    figure(2)
    plot(-K+1:K,pspec,'k','LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('k','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylabel('\eta(k,t)','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    ylim([-20 0])
    xlim([-256 256])
    % Plotting the paths in the four vortex case   
      
    figure(3)
    plot((0:count-1)*dt,xtrack(1,1:count),'k--',(0:count-1)*dt,xtrack(2,1:count),'b--',(0:count-1)*dt,xtrack(3,1:count),'k.',(0:count-1)*dt,xtrack(4,1:count),'b.','LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('t','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    
    figure(4)
    plot((0:count-1)*dt,ztrack(1,1:count),'k--',(0:count-1)*dt,ztrack(2,1:count),'b--',(0:count-1)*dt,ztrack(3,1:count),'k.',(0:count-1)*dt,ztrack(4,1:count),'b.','LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('t','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    
    % Plotting the paths in the two vortex case  
    %{  
    figure(3)
    plot((0:count-1)*dt,xtrack(1,1:count),'k--',(0:count-1)*dt,xtrack(2,1:count),'b--','LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('t','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    
    figure(4)
    plot((0:count-1)*dt,ztrack(1,1:count),'k--',(0:count-1)*dt,ztrack(2,1:count),'b--','LineWidth',2)
    set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
    xlabel('t','FontName','Helvetica','FontSize',30,'FontWeight','bold')
    %}
    % And if you want to make an animation of the surface uncomment this.  
      
    % Movie_Maker_1d(movie_plot,Xmesh,plot_count,'surf_response')
end