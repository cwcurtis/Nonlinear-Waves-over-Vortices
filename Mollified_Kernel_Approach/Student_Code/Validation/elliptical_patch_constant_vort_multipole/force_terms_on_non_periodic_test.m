function Kvec = force_terms_on_non_periodic_test(rval,gvals,xpos,zpos)

Nvorts = length(xpos);
Kx = zeros(Nvorts);
Kz = zeros(Nvorts);
Kvec = zeros(Nvorts,2);

for jj=1:Nvorts
    if(jj>1 && jj<Nvorts)
        dx = xpos(jj) - [xpos(1:jj-1);xpos(jj+1:end)];
        dzm = zpos(jj) - [zpos(1:jj-1);zpos(jj+1:end)];        
    elseif(jj==1)
        dx = xpos(1) - xpos(2:end);
        dzm = zpos(1) - zpos(2:end);        
    elseif(jj==Nvorts)
        dx = xpos(Nvorts) - xpos(1:Nvorts-1);
        dzm = zpos(Nvorts) - zpos(1:Nvorts-1);        
    end
     
    diff = dx.^2 + dzm.^2;
    er = exp(-diff/(2*rval^2));
    scv = (1-er).*(1+2*er)./diff;
    kerx = -dzm.*scv;
    kerz = dx.*scv;   
    if(jj>1 && jj<Nvorts)
        Kx(jj,1:jj-1) = kerx(1:jj-1);
        Kx(jj,jj+1:Nvorts) = kerx(jj:Nvorts-1);
        Kz(jj,1:jj-1) = kerz(1:jj-1);
        Kz(jj,jj+1:Nvorts) = kerz(jj:Nvorts-1);
    elseif(jj==1)
        Kx(1,2:Nvorts) = kerx;
        Kz(1,2:Nvorts) = kerz;
    elseif(jj==Nvorts)
        Kx(Nvorts,1:Nvorts-1) = kerx;
        Kz(Nvorts,1:Nvorts-1) = kerz;
    end
    
end

Kvec(:,1) = Kx*gvals;
Kvec(:,2) = Kz*gvals;