function rvec = linear_response(K,mu,gam,F,xv,zv,tf)

    KT = 2*K;
    
    Xmesh = linspace(-1,1,KT+1);
    Xmesh = Xmesh(1:KT)';    
    
    ep = mu*gam^2;
    
    omega = @(k) sqrt(ep*pi*k*tanh(pi*gam*k)/(gam*F^2));
    
    f = @(x,k) 2*sin(omega(k).*tf)./omega(k).*sinh(pi.*gam.*k.*zv)./gam.*exp(-pi.*gam.*k)...
               .*sin(pi.*k.*xv).*cos(pi.*k.*x);
    
    rvec = zeros(KT,1);
    
    for kk=1:30
       rvec = rvec + f(Xmesh,kk); 
    end
    
end

