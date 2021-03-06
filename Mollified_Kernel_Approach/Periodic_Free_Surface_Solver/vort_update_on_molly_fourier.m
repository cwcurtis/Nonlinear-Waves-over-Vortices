function [u,xpos,zpos] = vort_update_on_molly_fourier(Xmesh,gam,mu,ep,u,gvals,L1,no_dno_term,Nvorts,Ehdt,Edt,xpos,zpos,dt,sig,Mx,KTT)

k1 = dt*force_terms_on_molly_fourier(Xmesh,gam,mu,ep,u,gvals,L1,no_dno_term,sig,Mx,Nvorts);
uh1 = [Ehdt*(u(1:KTT)+k1(1:KTT)/2);xpos+k1(KTT+1:KTT+Nvorts)/2;zpos+k1(KTT+Nvorts+1:KTT+2*Nvorts)/2];

k2 = dt*force_terms_on_molly_fourier(Xmesh,gam,mu,ep,uh1,gvals,L1,no_dno_term,sig,Mx,Nvorts);
uh2 = [(Ehdt*u(1:KTT)+k2(1:KTT)/2);xpos+k2(KTT+1:KTT+Nvorts)/2;zpos+k2(KTT+Nvorts+1:KTT+2*Nvorts)/2];

k3 = dt*force_terms_on_molly_fourier(Xmesh,gam,mu,ep,uh2,gvals,L1,no_dno_term,sig,Mx,Nvorts);
uh3 = [(Edt*u(1:KTT)+Ehdt*k3(1:KTT));xpos+k3(KTT+1:KTT+Nvorts);zpos+k3(KTT+Nvorts+1:KTT+2*Nvorts)];

k4 = dt*force_terms_on_molly_fourier(Xmesh,gam,mu,ep,uh3,gvals,L1,no_dno_term,sig,Mx,Nvorts);

u(1:KTT) = Edt*u(1:KTT) + (Edt*k1(1:KTT) + 2*Ehdt*(k2(1:KTT) + k3(1:KTT)) + k4(1:KTT))/6;
u(KTT+1:KTT+2*Nvorts) = u(KTT+1:KTT+2*Nvorts) + ( k1(KTT+1:KTT+2*Nvorts) + 2*(k2(KTT+1:KTT+2*Nvorts)+k3(KTT+1:KTT+2*Nvorts)) + k4(KTT+1:KTT+2*Nvorts) )/6;

xpos = u(KTT+1:KTT+Nvorts);
zpos = u(KTT+Nvorts+1:KTT+2*Nvorts);

inds = xpos < -Mx;
if length(inds)>1
    xpos(inds) = xpos(inds) + 2*Mx;
end

inds = xpos > Mx;
if length(inds)>1
    xpos(inds) = xpos(inds) - 2*Mx;
end