function dnohot = dno_maker(eta,Q,phi0,Lop,gam,mu,Kmesh,Not)

    Dg = pi*gam*Kmesh;
    
    %del = mu;
    KT = length(Kmesh);
    K = KT/2;
    Kc = floor(2*K/3);
    Kuc = 2*K-Kc+1;
    Kc = Kc + 1;
    
    phis = zeros(KT,Not+1);    
    dnohot = zeros(KT,1);
    phis(:,1) = phi0;
    mup = 1;
    
    for jj=2:Not+1
  
        etap = ones(KT,1);
        Dgp = ones(KT,1);
        mup = mu*mup;
        
        for ll=1:jj-1            
                
                %etap = del*etap.*eta/ll;
                etap = etap.*eta/ll;
                Dgp = Dg.*Dgp;
                
                % De-alias while computing the dno.  Remove for speed-up.
                
                etap = fft(etap);
                etap(Kc:Kuc) = 0;
                etap = real(ifft(etap));
                        
                if(mod(ll,2)~=0)                
                    phis(:,jj) = phis(:,jj) -  gam^2*Dgp/(-1i*gam).*Lop.*fft(etap.*phis(:,jj-ll));            
                else            
                    phis(:,jj) = phis(:,jj) -  Dgp.*fft(etap.*phis(:,jj-ll));
                end                
        end    
        
        if(mod(jj,2)==0)                    
            phis(:,jj) = phis(:,jj) - Dgp/(-1i*gam).*fft(etap.*Q);
        else                      
            phis(:,jj) = phis(:,jj) - 1i*gam*Lop.*Dgp/(-1i*gam).*fft(etap.*Q);
        end        
        
        phihold = phis(:,jj);
        phihold(Kc:Kuc) = 0;
        phis(:,jj) = phihold;
        phis(:,jj) = real(ifft(phis(:,jj)));
        dnohot = dnohot + mup*phis(:,jj);
        
    end    
    
    %dnohot = sum(phis(:,2:Not+1),2);
        
end