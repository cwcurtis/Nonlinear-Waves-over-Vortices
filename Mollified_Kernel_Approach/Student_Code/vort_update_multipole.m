function [u,xpos,zpos] = vort_update_multipole(Xmesh,gam,mu,ep,F,u,gvals,L1,no_dno_term,Nvorts,Ntrunc,Ehdt,Edt,xpos,zpos,dt,KTT)

k1 = dt*force_terms_multipole(Xmesh,gam,mu,ep,F,u,gvals,L1,no_dno_term,Nvorts,Ntrunc);
uh1 = [Ehdt*(u(1:KTT)+k1(1:KTT)/2);xpos+k1(KTT+1:KTT+Nvorts)/2;zpos+k1(KTT+Nvorts+1:KTT+2*Nvorts)/2];

k2 = dt*force_terms_multipole(Xmesh,gam,mu,ep,F,uh1,gvals,L1,no_dno_term,Nvorts,Ntrunc);
uh2 = [(Ehdt*u(1:KTT)+k2(1:KTT)/2);xpos+k2(KTT+1:KTT+Nvorts)/2;zpos+k2(KTT+Nvorts+1:KTT+2*Nvorts)/2];

k3 = dt*force_terms_multipole(Xmesh,gam,mu,ep,F,uh2,gvals,L1,no_dno_term,Nvorts,Ntrunc);
uh3 = [(Edt*u(1:KTT)+Ehdt*k3(1:KTT));xpos+k3(KTT+1:KTT+Nvorts);zpos+k3(KTT+Nvorts+1:KTT+2*Nvorts)];

k4 = dt*force_terms_multipole(Xmesh,gam,mu,ep,F,uh3,gvals,L1,no_dno_term,Nvorts,Ntrunc);

u(1:KTT) = Edt*u(1:KTT) + (Edt*k1(1:KTT) + 2*Ehdt*(k2(1:KTT) + k3(1:KTT)) + k4(1:KTT))/6;
u(KTT+1:KTT+2*Nvorts) = u(KTT+1:KTT+2*Nvorts) + ( k1(KTT+1:KTT+2*Nvorts) + 2*(k2(KTT+1:KTT+2*Nvorts)+k3(KTT+1:KTT+2*Nvorts)) + k4(KTT+1:KTT+2*Nvorts) )/6;

xpos = u(KTT+1:KTT+Nvorts);
zpos = u(KTT+Nvorts+1:KTT+2*Nvorts);

% Relocate vortices that have crossed the boundary
for kk = 1:Nvorts
    while(xpos(kk) < -1)
        xpos(kk) = xpos(kk)+2;
    end

    while(xpos(kk) > 1)
        xpos(kk) = xpos(kk)-2;
    end
end