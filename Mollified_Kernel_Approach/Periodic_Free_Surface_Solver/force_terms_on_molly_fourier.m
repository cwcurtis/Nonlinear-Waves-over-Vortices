function nl = force_terms_on_molly_fourier(Xmesh,gam,mu,ep,u,gvals,L1,no_dno_term,sig,Mx,Nvorts)

KT = length(Xmesh);
K = KT/2;
%pM = pi/Mx;
p2M = pi/(2*Mx);

Kc = floor(2*K/3);
Kuc = 2*K-Kc+1;
Kc = Kc + 1;

Kmesh = 1/Mx*[0:K-1 0 -K+1:-1]';
Dx = pi*1i*Kmesh; 

eta = u(1:KT);
Q = u(KT+1:2*KT);
xpos = u(2*KT+1:2*KT+Nvorts);
zpos = u(2*KT+Nvorts+1:2*KT+2*Nvorts);

etax = real(ifft(Dx.*eta));
G0 = real(ifft(L1.*Q));        
eta = real(ifft(eta));
Q = real(ifft(Q));

dnonl = dno_maker(eta,Q,G0,L1,gam,mu,Kmesh,no_dno_term);
dno = G0 + dnonl; %update G_j

Pv = zeros(length(Xmesh),1);
Ev = zeros(length(Xmesh),1);

phix = zeros(length(Xmesh),1);
phiz = zeros(length(Xmesh),1);

%fphiz = @(x,z,zj) sin(pM*x).*sinh(pM*gam*z).*sinh(pM*gam*zj)./( ( cosh(pM*gam*(z-zj)) - cos(pM*x) ).*( cosh(pM*gam*(z+zj)) - cos(pM*x) ) );
%fphiza = @(x,z,zj) sin(pM*x).*(cosh(pM*gam*zj).*cosh(pM*gam*z) - cos(pM*x))./( ( cosh(pM*gam*(z-zj)) - cos(pM*x) ).*( cosh(pM*gam*(z+zj)) - cos(pM*x) ) );
%fphix = @(x,z,zj) sinh(pM*gam*zj).*(cosh(pM*gam*zj) - cosh(pM*gam*z).*cos(pM*x))./( (cosh(pM*gam*(z-zj)) - cos(pM*x)).*( cosh(pM*gam*(z+zj)) - cos(pM*x) ) );

%fpsix = @(x,z,zj) -.25*sin(pM*x).*( 1./( cosh(pi*gam*(z-zj)) - cos(pM*x) ) + 1./( cosh(pi*gam*(z+zj)) - cos(pM*x) ) );
%fpsiz = @(x,z,zj) -.25*( sinh(pM*gam*(z-zj))./( cosh(pM*gam*(z-zj)) - cos(pM*x) ) + sinh(pM*gam*(z+zj))./( cosh(pM*gam*(z+zj)) - cos(pM*x) ) );

%ftpsix = @(x,z,zj) -.25*sin(pM*x).*( 1./( cosh(pM*gam*(z-zj)) - cos(pM*x) ) - 1./( cosh(pM*gam*(z+zj)) - cos(pM*x) ) );
%ftpsiz = @(x,z,zj) -.25*( sinh(pM*gam*(z-zj))./( cosh(pM*gam*(z-zj)) - cos(pM*x) ) - sinh(pM*gam*(z+zj))./( cosh(pM*gam*(z+zj)) - cos(pM*x) ) );

pervel = @(x,z) 1/(4*Mx)*cot(p2M*(x+1i*gam*z));
fphiz = @(x,z,zj) real(pervel(x,z-zj)-pervel(x,z+zj));
fphiza = @(x,z,zj) real(pervel(x,z-zj)+pervel(x,z+zj));
fphix = @(x,z,zj) imag(pervel(x,z-zj)-pervel(x,z+zj));

fpsix = @(x,z,zj) -real(pervel(x,z-zj)+pervel(x,z+zj));
fpsiz = @(x,z,zj) imag(pervel(x,z-zj)+pervel(x,z+zj));

ftpsix = @(x,z,zj) -real(pervel(x,z-zj)-pervel(x,z+zj));
ftpsiz = @(x,z,zj) imag(pervel(x,z-zj)-pervel(x,z+zj));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if Nvorts < 2024
%    veloc = direct_solver_periodic(ep,xpos,gam*zpos,gvals,Mx);
%else
    pval = 10;
    tree_vals = multi_pole_kernel_build(xpos,gam*zpos,gvals,pval,Mx,Nvorts);
    tree_vals = multi_pole_list_maker(tree_vals,pval,Nvorts);
    veloc = multi_pole_kernel_quick(xpos,gam*zpos,gvals,ep,pval,Mx,tree_vals);    
%end

xdot = 1/(4*Mx)*veloc(:,1);
zdot = 1/(4*Mx)*veloc(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj=1:Nvorts   
    
    psix = fpsix(Xmesh-xpos(jj),1+mu*eta,zpos(jj));
    psiz = fpsiz(Xmesh-xpos(jj),1+mu*eta,zpos(jj));
    
    tpsix = ftpsix(Xmesh-xpos(jj),1+mu*eta,zpos(jj));
    tpsiz = ftpsiz(Xmesh-xpos(jj),1+mu*eta,zpos(jj));
        
    bp1 = 2*Mx*sum(gam*dno.*psix + Q.*psiz)/KT;
    bp2 = 2*Mx*sum(gam*dno.*tpsiz - Q.*tpsix)/KT; 
    
    xdot(jj) = real(mu*xdot(jj) - mu*bp1);
    zdot(jj) = real(mu/gam*zdot(jj) - mu/gam*bp2);    
end

for jj=1:length(Xmesh)
    xvec = gvals.*( fphix(Xmesh(jj)-xpos,1 + mu*eta(jj),zpos) );
    zvec = gvals.*( fphiz(Xmesh(jj)-xpos,1 + mu*eta(jj),zpos) );
    zveca = gvals.*( fphiza(Xmesh(jj)-xpos,1 + mu*eta(jj),zpos) );
    
    phix(jj) = sum( xvec );
    phiz(jj) = sum( zvec );
    Pv(jj) = phiz(jj) - mu*gam*etax(jj).*phix(jj);
    Ev(jj) = sum((xdot.*xvec + gam*zdot.*zveca))-mu*(phix(jj)^2 + phiz(jj)^2)/2;    
end
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