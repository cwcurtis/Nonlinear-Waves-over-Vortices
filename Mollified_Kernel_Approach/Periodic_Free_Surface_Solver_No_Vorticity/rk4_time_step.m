function u = rk4_time_step(Xmesh,Mx,gam,mu,u,L1,no_dno_term,Ehdt,Edt,dt)

k1 = dt*force_terms(Xmesh,Mx,gam,mu,u,L1,no_dno_term);
uh1 = Ehdt*(u+k1/2);

k2 = dt*force_terms(Xmesh,Mx,gam,mu,uh1,L1,no_dno_term);
uh2 = Ehdt*u+k2/2;

k3 = dt*force_terms(Xmesh,Mx,gam,mu,uh2,L1,no_dno_term);
uh3 = Edt*u+Ehdt*k3;

k4 = dt*force_terms(Xmesh,Mx,gam,mu,uh3,L1,no_dno_term);

u = Edt*u + (Edt*k1 + 2*Ehdt*(k2 + k3) + k4)/6;

