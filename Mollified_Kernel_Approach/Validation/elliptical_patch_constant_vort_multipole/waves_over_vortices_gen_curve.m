function [tvals,evec] = waves_over_vortices_gen_curve(Nx,mu,gam,omega,tf,samp)
    close all
    
    % Choose time step and find inverse of linear part of semi-implicit
    % time stepping scheme.
    Rv = 1/10;
    zoff = .5;
    axs = 1;
    av = Rv*axs;
    bv = Rv;
    
    F = pi*omega*av*bv/gam;
    pmesh = linspace(0,2*pi,2*Nx);
     
    cp = cos(pmesh);
    sp = sin(pmesh);    
    %tdot = mu*omega*av*bv/(av + gam*bv)^2;
    
    [xpos,zpos,gvals,rval,Nvorts] = initializer(Nx,omega,gam,av,bv,1);
    
    dt = .05;
    nmax = round(tf/dt);
    
    zpos = zoff + zpos/gam;
    
    %inter = 5;
    plot_count = 1;
    no_of_evals = round(nmax/samp);
    Vcnt = zeros(no_of_evals+1,1);
    Vcnt(1) = Nvorts;
    xtrack = xpos;
    ztrack = zpos;
    gtrack = gvals;
    
    u = [xpos;zpos]; %velocity vector field
    
    evec = zeros(no_of_evals,1);
    tvals = zeros(no_of_evals,1);
    cnt = 1;
    
    %Compute reference vorticity profile
    %{
    chi = @(r) 1./pi.*(2.*exp(-(r.^2))-.5*exp(-(r/sqrt(2)).^2));
    
    imat = zeros(length(zvals),length(xvals));
    
    for jj=1:length(xvals)
       for kk=1:length(zvals)
          diff = sqrt((xvals(jj)-xpos).^2 + gam.^2*(zvals(kk)-zpos).^2)/rval;          
          imat(kk,jj) = gam*sum(gvals.*chi(diff)/(rval^2));
       end
    end   
    disp(sqrt(sum(sum(abs(imat-rmat).^2))/sum(sum(rmat.^2))))
    dmat = log10(abs(imat-rmat));
    figure(1)
    surf(xvals,zvals,dmat,'LineStyle','none')
    figure(2)
    plot(xvals,dmat(Nvx/2+1,:),'LineWidth',2)
    pause
    %}
    for jj=1:nmax      
        % Now update the vortex positions         
        u = vort_update_multipole(mu,gam,rval,u,gvals,Nvorts,dt);           
            
        if(mod(jj,samp)==0)
            plot_count = plot_count + 1;
            
            xtrack = [xtrack;xpos];
            ztrack = [ztrack;zpos];
            gtrack = [gtrack;gvals];
            Vcnt(plot_count) = Nvorts;            
                    
            [xpos,zpos,gvals] = recircer(gvals,u(1:Nvorts),u(Nvorts+1:2*Nvorts),Nx);
            Nvorts = length(gvals);
            disp('Number of Vortices is')
            disp(Nvorts)
            u = zeros(2*Nvorts,1);
            u(1:Nvorts) = xpos;
            u(Nvorts+1:2*Nvorts) = zpos;             
            
            evec(cnt) = interp_error(xpos,zpos,gvals,rval,gam,omega,av,bv,zoff);
            tvals(cnt) = dt*cnt*samp;
            cnt = cnt + 1;
        end
    end     
       
    %figure(2)
    %plot(tvals,evec,'k-','LineWidth',2)    
    %h = set(gca,'FontSize',30);
    %set(h,'Interpreter','LaTeX')
    %xlabel('$t$','Interpreter','LaTeX','FontSize',30)
    %ylabel('$\mathcal{E}(t)$','Interpreter','LaTeX','FontSize',30)
end