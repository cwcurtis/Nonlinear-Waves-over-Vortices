import numpy as np

def waves_over_vortices_gen_curve(Nx,mu,gam,omega,tf):
    Rv = 1./25.
    zoff = .25
    axs = 2.
    av = Rv*axs
    bv = Rv
    
    F = pi*omega*av*bv/gam
    Gamma = F
    pmesh = np.linspace(0,2*pi,2*Nx)
     
    cp = np.cos(pmesh)
    sp = np.sin(pmesh)    
    
    [xpos,zpos,gvals,rval,Nvorts] = initializer(Nx,omega,gam,av,bv)
    
    dt = .1
    nmax = int(np.round(tf/dt))
    
    xtrack = xpos
    ztrack = zoff + zpos/gam
    gtrack = gvals
    
    inter = 2
    plot_count = 0
    no_of_evals = round(nmax/inter)
    times = np.zeros(no_of_evals)
    errors = np.zeros(no_of_evals)
    Vcnt = np.zeros(no_of_evals+1)
    Vcnt[0] = Nvorts
    
    u = np.concatenate((xpos,zpos),axis=0) 
    
    for jj in xrange(1,nmax+1):
      
        u = vort_update_multipole(mu,1,rval,u,gvals,Nvorts,dt)   
         
        if(np.mod(jj,inter)==0):
            times(plot_count) = (jj-1)*dt
            xpos = u(1:Nvorts)
            zpos = zoff + u(Nvorts+1:2*Nvorts)/gam
            
            tv = mu/gam*omega*av*bv/(av+bv)^2*(jj-1)*dt
            ca = np.cos(tv)
            sa = np.sin(tv)
            xxs = cp*av
            yxs = sp*bv
            xellip = xxs*ca - yxs*sa
            yellip = (xxs*sa + yxs*ca)/gam + zoff
            
            #error = interp_error(xpos,zpos,gvals,rval,gam,omega,av,bv,tv,zoff);
            #errors(plot_count) = error;
                        
            xtrack = np.concatenate((xtrack,xpos),axis=0)
            ztrack = np.concatenate((ztrack,zpos),axis=0)
            gtrack = np.concatenate((gtrack,gvals),axis=0)
            plot_count += 1
            
            Vcnt[plot_count] = Nvorts   
            
           
        #if(np.mod(jj,2*inter)==0):
        #    [xpud,zpud,gvud] = recircer(gvals,xpos,gam*(zpos-zoff),Nx);
        #    Nvorts = length(gvud);
        #    gvals = gvud;
        #    u = zeros(2*Nvorts,1);
        #    u(1:Nvorts) = xpud;
        #    u(Nvorts+1:2*Nvorts) = zpud;            
        
            
    
    