function nl = force_terms_on_molly_fourier(Xmesh,Mx,gam,mu,u,L1,no_dno_term)

KT = length(Xmesh);
K = KT/2;

Kc = floor(2*K/3);
Kuc = 2*K-Kc+1;
Kc = Kc + 1;

Kmesh = [0:K-1 0 -K+1:-1]';
Dx = pi*1i*Kmesh/Mx; 

eta = u(1:KT);
Q = u(KT+1:2*KT);

etax = real(ifft(Dx.*eta));
G0 = real(ifft(L1.*Q));        
eta = real(ifft(eta));
Q = real(ifft(Q));

dnonl = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh,no_dno_term);
dno = G0 + dnonl; %update G_j

gm = mu*gam;

% De-alias the terms we need to find the nonlinearity associated with the
% Bernoulli equation.

term1 = etax.*Q;
term1 = fft(term1);
term1(Kc:Kuc) = 0;
term1 = real(ifft(term1));

nl2 = -mu./(1+(gm*etax).^2).*( .5*Q.^2 - gm*gam*dno.*term1 - .5*(gam*dno).^2 );

nl1 = fft(dnonl);
nl2 = Dx.*fft(nl2);

% De-alias Bernoulli equation nonlinearity on our way out of this
% subroutine.

nl1(Kc:Kuc) = 0;
nl2(Kc:Kuc) = 0;

nl = [nl1;nl2];