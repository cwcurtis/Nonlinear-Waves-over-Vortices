function waves_over_vortices_gen_curve(Nx,mu,gam,omega,tf,npow)
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
    
    [xpos,zpos,gvals,rval,Nvorts] = initializer(Nx,omega,gam,av,bv,npow);
    
    dt = .05;
    nmax = round(tf/dt);
    
    zpos = zoff + zpos/gam;
    
    inter = 5;
    plot_count = 1;
    no_of_evals = round(nmax/inter);
    Vcnt = zeros(no_of_evals+1,1);
    Vcnt(1) = Nvorts;
    xtrack = xpos;
    ztrack = zpos;
    gtrack = gvals;
    
    u = [xpos;zpos]; %velocity vector field
    
    % Make folder
    S = make_folder(Nx/2,Nx,mu,gam,F,tf,av/bv);
    clf
    n_bdry = 0;     % Number of points in cicular boundary.
    markersize = 10;
    evec = zeros(nmax/inter,1);
    tvals = zeros(nmax/inter,1);
    cnt = 1;
    
    %Compute reference vorticity profile
    Nvx = 200;
    xvals = linspace(-.2,.2,Nvx)';
    zvals = linspace(.1,.9,Nvx)';
    rmat = zeros(length(zvals),length(xvals));
    cfun = @(x,z) (x/av).^2 + gam^2*(z/bv).^2;
    ofun = @(x,z) omega*(1-cfun(x,z)).^3;
    for ll = 1:Nvx
       inds = cfun(xvals(ll),zvals-zoff) <= 1;
       if length(inds)>1
           rmat(ll,inds) = ofun(xvals(ll),zvals(inds)-zoff); 
       end
    end
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
            
        if(mod(jj,inter)==0)
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
            
            evec(cnt) = interp_error(xpos,zpos,gvals,rval,gam,omega,av,bv,zoff,rmat);
            tvals(cnt) = dt*cnt*inter;
            cnt = cnt + 1;
        end
    end     
     
    figure(1)
    gif_my_gif(xtrack,ztrack,gtrack,n_bdry,Vcnt,plot_count,S,markersize);    
    
    figure(2)
    plot(tvals,evec,'k-','LineWidth',2)    
    
end