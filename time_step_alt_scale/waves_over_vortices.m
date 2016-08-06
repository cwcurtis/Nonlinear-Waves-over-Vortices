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
    
    dt = 1e-3;
    nmax = round(tf/dt);
    
    L1 = -1i*tanh(pi*gam*Kmesh)/gam;
    L2 = -1i*pi*Kmesh/F^2;
    
    Linv11 = (1 - 9/16*dt^2*L1.*L2).^(-1);
    Linv22 = Linv11;
    Linv12 = 3/4*dt*L1.*Linv11;
    Linv21 = 3/4*dt*L2.*Linv11;
    
    % Here we put in a flat surface and zero background velocity potential
    %{
    etan = zeros(KT,1);
    qn = zeros(KT,1);    
    %}
    % Here we put in the cnoidal solution for KdV in the appropriate
    % coordinates.  
    
    modu = .2;    
    kap = ellipke(modu);
    uvals = kap*Xmesh;
    [~,cn,~] = ellipj(uvals,modu);
    
    ceff = F + 2/3*mu*F^2*kap^2*(2*modu^2-1);
    
    qn = fft(4/(3*F)*modu^2*kap^2*cn.^2);
    
    etan = ceff*qn;
    
    etan1 = etan;
    qn1 = qn;
    
    xtrack = zeros(Nvorts,nmax+1);
    ztrack = zeros(Nvorts,nmax+1);
    
    xtrack(:,1) = xpos;
    ztrack(:,1) = zpos;
    
    inter = 10; % Number of time steps between frames of movie
    no_of_plots = round(nmax/inter); % Number of plots in movie
    movie_plot = zeros(no_of_plots,KT); % Movie storage
    plot_count = 1;
    count = 1;
    
    no_dno_term = 2;
    
    for jj=1:nmax
        
        if(jj==1)
            
            etanr = real(ifft(etan));
            Q = real(ifft(qn));
            G0 = real(ifft(L1.*qn));
            etax = real(ifft(Dx.*etan));
            
            dnohot = dno_maker(etanr,Q,G0,L1,gam,mu,Kmesh,no_dno_term);
            [Pv,nl2n] = force_terms(Xmesh,gam,mu,etanr,etax,Q,xpos,zpos,gvals,G0+dnohot);        
            nl1n = fft(dnohot) + Pv/gam;
                
            nl1n1 = nl1n;
            nl1n2 = nl1n1;
            nl1n3 = nl1n2;
          
            nl2n1 = nl2n;
            nl2n2 = nl2n1;
            nl2n3 = nl2n2;
                        
        else
        
            etax = real(ifft(Dx.*etan));
            [Pv,nl2n] = force_terms(Xmesh,gam,mu,etanr,etax,Q,xpos,zpos,gvals,G0+dnohot);        
            nl1n = fft(dnohot) + Pv/gam;
            
        end
        
        nlvec1 = dt*(55/24*nl1n - 59/24*nl1n1 + 37/24*nl1n2 - 3/8*nl1n3);
        nlvec2 = dt*(55/24*nl2n - 59/24*nl2n1 + 37/24*nl2n2 - 3/8*nl2n3);
        
        etanp1 = -etan1/3 + Linv11.*(etan + etan1/3 + nlvec1) + Linv12.*(qn + qn1/3 + nlvec2);
        qnp1 = -qn1/3 + Linv21.*(etan + etan1/3 + nlvec1) + Linv22.*(qn + qn1/3 + nlvec2);
       
        % Interpolate to find the half-time surface and potential profiles
        etanph = -1/8*etan1 + 3/4*etan + 3/8*etanp1;
        qnph = -1/8*qn1 + 3/4*qn + 3/8*qnp1;
        
        % De-Alias
        etanph(Kc:Kuc) = 0;
        qnph(Kc:Kuc) = 0;
        
        etanrh = real(ifft(etanph));
        Qh = real(ifft(qnph));
        G0h = real(ifft(L1.*qnph));
        
        dnohoth = dno_maker(etanrh,Qh,G0h,L1,gam,mu,Kmesh,no_dno_term);
        
        etanp1(Kc:Kuc) = 0;
        qnp1(Kc:Kuc) = 0;
                    
        etanrp1 = real(ifft(etanp1));
        Qp1 = real(ifft(qnp1));
        G0p1 = real(ifft(L1.*qnp1));
        
        dnohotp1 = dno_maker(etanrp1,Qp1,G0p1,L1,gam,mu,Kmesh,no_dno_term);
        % Now update the vortex positions
        
        k1 = dt*vort_update(Xmesh,gam,mu,etanr,Q,xpos,zpos,gvals,G0+dnohot);
        k2 = dt*vort_update(Xmesh,gam,mu,etanrh,Qh,xpos+k1(1:Nvorts)/2,zpos+k1(Nvorts+1:2*Nvorts)/2,gvals,G0h+dnohoth);
        k3 = dt*vort_update(Xmesh,gam,mu,etanrh,Qh,xpos+k2(1:Nvorts)/2,zpos+k2(Nvorts+1:2*Nvorts)/2,gvals,G0h+dnohoth);
        k4 = dt*vort_update(Xmesh,gam,mu,etanrp1,Qp1,xpos+k3(1:Nvorts),zpos+k3(Nvorts+1:2*Nvorts),gvals,G0p1+dnohotp1);
        
        hvec = [xpos;zpos] + (k1 + 2*k2 + 2*k3 + k4)/6;
        xpos = hvec(1:Nvorts);
        zpos = hvec(Nvorts+1:2*Nvorts);
        
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
        
            etan1 = etan;
            qn1 = qn;
            etan = etanp1;
            qn = qnp1;
            
            G0 = G0p1;
            dnohot = dnohotp1;
            etanr = etanrp1;
            Q = Qp1;
            
            nl1n3 = nl1n2;
            nl1n2 = nl1n1;
            nl1n1 = nl1n;
          
            nl2n3 = nl2n2;
            nl2n2 = nl2n1;
            nl2n1 = nl2n;
            
            count = count + 1;
        end
        
        % Add to the movie plot matrix for the purpose of making animations
        
        if(mod(jj,inter)==0)
            movie_plot(plot_count,:) = real(ifft(etan));
            plot_count = plot_count + 1;
        end
        
    end
    
    etan(Kc:Kuc) = 0;
    qn(Kc:Kuc) = 0;
    nfin = real(ifft(etan));
    pspec = log10(fftshift(abs(etan))/KT);
    
    %qfin = real(ifft(qn));
    
    % If you want to compare against the analytic formula for the linear
    % motion
    %{
    linrep = linear_response(K,mu,gam,F,xpos(2),zpos(2),tf);
    figure(1)
    plot(Xmesh,nfin,'r--',Xmesh,linrep,'k.')
    %}
    
    % Plot the surface profile
    
    figure(1)
    plot(Xmesh,nfin)
       
    % Plot the power spectrum
    
    figure(2)
    plot(-K+1:K,pspec)
    
    % Plotting the paths in the four vortex case   
      
    
    figure(3)
    plot((0:count-1)*dt,xtrack(1,1:count),'k--',(0:count-1)*dt,xtrack(2,1:count),'b--',(0:count-1)*dt,xtrack(3,1:count),'k.',(0:count-1)*dt,xtrack(4,1:count),'b.')
    figure(4)
    plot((0:count-1)*dt,ztrack(1,1:count),'k--',(0:count-1)*dt,ztrack(2,1:count),'b--',(0:count-1)*dt,ztrack(3,1:count),'k.',(0:count-1)*dt,ztrack(4,1:count),'b.')
    
    
    % Plotting the paths in the two vortex case  
    %{  
    figure(3)
    plot((0:count-1)*dt,xtrack(1,1:count),'k--',(0:count-1)*dt,xtrack(2,1:count),'b--')
    figure(4)
    plot((0:count-1)*dt,ztrack(1,1:count),'k--',(0:count-1)*dt,ztrack(2,1:count),'b--')
    %}
    % And if you want to make an animation of the surface uncomment this.  
      
    %Movie_Maker_1d(movie_plot,Xmesh,plot_count,'surf_response')
end

