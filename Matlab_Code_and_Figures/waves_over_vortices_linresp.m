function tl = waves_over_vortices_linresp(K,mu,gam,F,tf)

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
    
    xpos = [-2*mu*gam;2*mu*gam];
    
    zpos1 = .25;
    zpos = [zpos1;zpos1];
       
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.  
    
    dt = 1e-2;
    nmax = round(tf/dt);
    
    L1 = -1i*tanh(pi.*gam.*Kmesh)./gam;
    L2 = -1i*pi*Kmesh;
    
    omf = @(k) 1i*sqrt(pi.*k.*tanh(pi.*gam.*k)./gam);
    vmatf = @(k) sqrt(pi.*gam.*k./tanh(pi.*gam.*k));
    vmatv = sign(Kmesh).*vmatf(Kmesh);
    vmatv(1) = 1;
    vmatv(K+1) = 1;
    
    Vmat = [eye(KT) eye(KT);-diag(vmatv) diag(vmatv)];
    Vmati = 1/2*[eye(KT) -diag(vmatv.^(-1)); eye(KT) diag(vmatv.^(-1))];
    
    Edt = Vmat*diag([exp(dt*omf(Kmesh));exp(-dt*omf(Kmesh))])*Vmati;
    Ehdt = Vmat*diag([exp(dt/2*omf(Kmesh));exp(-dt/2*omf(Kmesh))])*Vmati;
    
    % Here we put in a flat surface and zero background velocity potential    
    
    % Build the quiescent initial velocity potential 
    eta = zeros(KT,1);
    phix = init_cond(Xmesh,gam,xpos,zpos,gvals);
    Q = -F*fft(phix);
    Q(Kc:Kuc) = 0;
    Q0 = Q;
    
    xtrack = zeros(Nvorts,nmax+1);
    ztrack = zeros(Nvorts,nmax+1);
    
    xtrack(:,1) = xpos;
    ztrack(:,1) = zpos;
    
    no_dno_term = 20;
    
    u = [eta;Q;xpos;zpos];
    tl = tf;
    
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
        end
        
        linrep = linear_response(K,Q0,gam,F,xpos(2),zpos(2),jj*dt);
        eta = real(ifft(u(1:KT)));
        if(norm(linrep-eta)/norm(eta)>mu)
            tl = jj*dt;
            break
        end        
    end 
end