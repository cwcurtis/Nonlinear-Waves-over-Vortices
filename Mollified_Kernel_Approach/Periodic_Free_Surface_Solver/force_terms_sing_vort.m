function nl = force_terms_sing_vort(Xmesh,gam,mu,u,gval,L1,no_dno_term,sig,Mx)

KT = length(Xmesh);
K = KT/2;
p2M = pi/(2*Mx);

Kc = floor(2*K/3);
Kuc = 2*K-Kc+1;
Kc = Kc + 1;

Kmesh = 1/Mx*[0:K-1 0 -K+1:-1]';
Dx = pi*1i*Kmesh; 

eta = u(1:KT);
Q = u(KT+1:2*KT);
xpos = u(2*KT+1);
zpos = u(2*KT+2);

etax = real(ifft(Dx.*eta));
G0 = real(ifft(L1.*Q));        
eta = real(ifft(eta));
Q = real(ifft(Q));

dnonl = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh,no_dno_term);
dno = G0 + dnonl; %update G_j

pervel = @(x,z) 1/(4*Mx)*cot(p2M*(x+1i*gam*z));
fphiz = @(x,z,zj) real(pervel(x,z-zj)-pervel(x,z+zj));
fphiza = @(x,z,zj) real(pervel(x,z-zj)+pervel(x,z+zj));
fphix = @(x,z,zj) imag(pervel(x,z-zj)-pervel(x,z+zj));

fpsix = @(x,z,zj) -real(pervel(x,z-zj)+pervel(x,z+zj));
fpsiz = @(x,z,zj) imag(pervel(x,z-zj)+pervel(x,z+zj));

ftpsix = @(x,z,zj) -real(pervel(x,z-zj)-pervel(x,z+zj));
ftpsiz = @(x,z,zj) imag(pervel(x,z-zj)-pervel(x,z+zj));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
veloc = direct_solver_sing(gam*zpos,gval,Mx);

xdot = 1/(4*Mx)*veloc(:,1);
zdot = 1/(4*Mx)*veloc(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psix = fpsix(Xmesh-xpos,1+mu*eta,zpos);
psiz = fpsiz(Xmesh-xpos,1+mu*eta,zpos);
    
tpsix = ftpsix(Xmesh-xpos,1+mu*eta,zpos);
tpsiz = ftpsiz(Xmesh-xpos,1+mu*eta,zpos);
        
bp1 = 2*Mx*sum(gam*dno.*psix + Q.*psiz)/KT;
bp2 = 2*Mx*sum(gam*dno.*tpsiz - Q.*tpsix)/KT; 
    
xdot = real(mu*xdot - mu*bp1);
zdot = real(mu/gam*zdot - mu/gam*bp2);    

phix = gval*( fphix(Xmesh-xpos,1 + mu*eta,zpos) );
phiz = gval*( fphiz(Xmesh-xpos,1 + mu*eta,zpos) );
zveca = gval*( fphiza(Xmesh-xpos,1 + mu*eta,zpos) );
    
Pv = phiz - mu*gam*etax.*phix;
Ev = xdot*phix + gam*zdot*zveca - mu*(phix.^2 + phiz.^2)/2;    

gm = mu*gam;

% De-alias the terms we need to find the nonlinearity associated with the
% Bernoulli equation.

Pv = fft(Pv);
Pv(Kc:Kuc) = 0;
Pv = real(ifft(Pv));

phix = fft(phix);
phix(Kc:Kuc) = 0;
phix = real(ifft(phix));

phiz = fft(phiz);
phiz(Kc:Kuc) = 0;
phiz = real(ifft(phiz));

term1 = etax.*Q;
term1 = fft(term1);
term1(Kc:Kuc) = 0;
term1 = real(ifft(term1));

stterm = real(ifft(Dx.*fft(etax./sqrt(1+(gm*etax).^2))));

nl2 = (Ev - mu./(1+(gm*etax).^2).*( .5*Q.^2 - mu*gam*(gam*dno+Pv).*term1 +...
                                     phix.*(Q - gam*gm*etax.*dno) - .5*gam*dno.*(2*Pv + gam*dno) + phiz.*(gm*term1+gam*dno) ) + sig*stterm);

nl1 = fft(dnonl + Pv/gam);
nl2 = Dx.*fft(nl2);

% De-alias Bernoulli equation nonlinearity on our way out of this
% subroutine.

nl1(Kc:Kuc) = 0;
nl2(Kc:Kuc) = 0;

nl = [nl1;nl2;xdot;zdot];