function rvec = linear_response(K,qn0,gam,F,xv,zv,tf)

    KT = 2*K;
    
    Xmesh = linspace(-1,1,KT+1);
    Xmesh = Xmesh(1:KT)';    
    Kmesh = [0:K -K+1:-1]';
    
    omega = @(k) sqrt(pi*k.*tanh(pi*gam*k)/gam);
    
    f = @(x,k) 2*sin(omega(k).*tf)./omega(k).*sinh(pi.*gam.*k.*zv)./gam.*exp(-pi.*gam.*k)...
               .*sin(pi.*k.*xv).*cos(pi.*k.*x);
    
    rvec = zeros(KT,1);
    
    for kk=1:30
       rvec = rvec + f(Xmesh,kk); 
    end
    
    icvec = omega(Kmesh)./(1i*pi*Kmesh).*sin(omega(Kmesh)*tf);
    icvec(1) = 0;
    %disp(size(icvec))
    %disp(size(qn))
    rvec = F*rvec + real(ifft(qn0.*icvec));
    
end

