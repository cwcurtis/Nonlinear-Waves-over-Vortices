function [Edt,Ehdt] = get_stuff(Kmesh,gam,K,KT,dt)

omf = @(k) 1i*sqrt(pi.*k.*tanh(pi.*gam.*k)./gam);
disp('Stability parameter for RK4:')
disp(dt*max(abs(omf(Kmesh))))
vmatf = @(k) sqrt(pi.*gam.*k./tanh(pi.*gam.*k));
vmatv = sign(Kmesh).*vmatf(Kmesh);
vmatv(1) = 1;
vmatv(K+1) = 1;

Vmat = [eye(KT) eye(KT);-diag(vmatv) diag(vmatv)];
Vmati = 1/2*[eye(KT) -diag(vmatv.^(-1)); eye(KT) diag(vmatv.^(-1))];

Edt = Vmat*diag([exp(dt*omf(Kmesh));exp(-dt*omf(Kmesh))])*Vmati;
Ehdt = Vmat*diag([exp(dt/2*omf(Kmesh));exp(-dt/2*omf(Kmesh))])*Vmati;