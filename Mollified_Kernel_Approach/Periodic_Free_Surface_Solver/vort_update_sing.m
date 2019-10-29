function [u,xpos,zpos] = vort_update_sing(Xmesh,gam,mu,u,gval,L1,no_dno_term,Ehdt,Edt,xpos,zpos,dt,sig,Mx,KTT)

k1 = dt*force_terms_sing_vort(Xmesh,gam,mu,u,gval,L1,no_dno_term,sig,Mx);
uh1 = [Ehdt*(u(1:KTT)+k1(1:KTT)/2);xpos+k1(KTT+1)/2;zpos+k1(KTT+2)/2];

k2 = dt*force_terms_sing_vort(Xmesh,gam,mu,uh1,gval,L1,no_dno_term,sig,Mx);
uh2 = [(Ehdt*u(1:KTT)+k2(1:KTT)/2);xpos+k2(KTT+1)/2;zpos+k2(KTT+2)/2];

k3 = dt*force_terms_sing_vort(Xmesh,gam,mu,uh2,gval,L1,no_dno_term,sig,Mx);
uh3 = [(Edt*u(1:KTT)+Ehdt*k3(1:KTT));xpos+k3(KTT+1);zpos+k3(KTT+2)];

k4 = dt*force_terms_sing_vort(Xmesh,gam,mu,uh3,gval,L1,no_dno_term,sig,Mx);

u(1:KTT) = Edt*u(1:KTT) + (Edt*k1(1:KTT) + 2*Ehdt*(k2(1:KTT) + k3(1:KTT)) + k4(1:KTT))/6;
u(KTT+1:KTT+2) = u(KTT+1:KTT+2) + ( k1(KTT+1:KTT+2) + 2*(k2(KTT+1:KTT+2)+k3(KTT+1:KTT+2)) + k4(KTT+1:KTT+2) )/6;

xpos = u(KTT+1);
zpos = u(KTT+2);

inds = xpos < -Mx;
if length(inds)>1
    xpos(inds) = xpos(inds) + 2*Mx;
end

inds = xpos > Mx;
if length(inds)>1
    xpos(inds) = xpos(inds) - 2*Mx;
end