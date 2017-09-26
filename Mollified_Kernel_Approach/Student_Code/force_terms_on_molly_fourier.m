function nl = force_terms_on_molly_fourier(Xmesh,gam,mu,ep,F,u,gvals,L1,no_dno_term,Nvorts,Ntrunc)

KT = length(Xmesh);
K = KT/2;

Kc = floor(2*K/3);
Kuc = 2*K-Kc+1;
Kc = Kc + 1;

Kmesh = [0:K-1 0 -K+1:-1]';
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

xdot = zeros(Nvorts,1);
zdot = zeros(Nvorts,1);

Pv = zeros(length(Xmesh),1);
Ev = zeros(length(Xmesh),1);

phix = zeros(length(Xmesh),1);
phiz = zeros(length(Xmesh),1);

fphiz = @(x,z,zj) sin(pi*x).*sinh(pi*gam*z).*sinh(pi*gam*zj)./( ( cosh(pi*gam*(z-zj)) - cos(pi*x) ).*( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );
fphiza = @(x,z,zj) sin(pi*x).*(cosh(pi*gam*zj).*cosh(pi*gam*z) - cos(pi*x))./( ( cosh(pi*gam*(z-zj)) - cos(pi*x) ).*( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );
fphix = @(x,z,zj) sinh(pi*gam*zj).*(cosh(pi*gam*zj) - cosh(pi*gam*z).*cos(pi*x))./( (cosh(pi*gam*(z-zj)) - cos(pi*x)).*( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );

fpsix = @(x,z,zj) -.25*sin(pi*x).*( 1./( cosh(pi*gam*(z-zj)) - cos(pi*x) ) + 1./( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );
fpsiz = @(x,z,zj) -.25*( sinh(pi*gam*(z-zj))./( cosh(pi*gam*(z-zj)) - cos(pi*x) ) + sinh(pi*gam*(z+zj))./( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );

ftpsix = @(x,z,zj) -.25*sin(pi*x).*( 1./( cosh(pi*gam*(z-zj)) - cos(pi*x) ) - 1./( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );
ftpsiz = @(x,z,zj) -.25*( sinh(pi*gam*(z-zj))./( cosh(pi*gam*(z-zj)) - cos(pi*x) ) - sinh(pi*gam*(z+zj))./( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );

Kx = zeros(Nvorts);
Kz = zeros(Nvorts);

parfor jj=1:Nvorts
    dx = xpos(jj) - xpos;
    dzm = zpos(jj) - zpos;
    dzp = zpos(jj) + zpos;
    for ll=1:Nvorts                
        % compute the mollified part of the kernel
        [msumxm,msumzm] = kernel_mol(dx(ll),dzm(ll),gam,ep,Ntrunc);
        fac = cosh(gam*pi*dzp(ll)) - cos(pi*dx(ll));
        msumxp = -sinh(gam*pi*dzp(ll))/fac;
        msumzp = sin(pi*dx(ll))/fac;
        
        Kx(ll,jj) = msumxm-msumxp;
        Kz(ll,jj) = msumzm-msumzp;
    end
end


for jj=1:Nvorts
   
    xdot(jj) = mu*F*sum(gvals.*Kx(:,jj));
    zdot(jj) = mu*F/gam*sum(gvals.*Kz(:,jj));
    
    psix = fpsix(Xmesh-xpos(jj),1+mu*eta,zpos(jj));
    psiz = fpsiz(Xmesh-xpos(jj),1+mu*eta,zpos(jj));
    
    tpsix = ftpsix(Xmesh-xpos(jj),1+mu*eta,zpos(jj));
    tpsiz = ftpsiz(Xmesh-xpos(jj),1+mu*eta,zpos(jj));
        
    bp1 = 2*sum(gam*dno.*psix + Q.*psiz)/KT;
    bp2 = 2*sum(gam*dno.*tpsiz - Q.*tpsix)/KT; 
    
    xdot(jj) = real(xdot(jj) - mu*bp1);
    zdot(jj) = real(zdot(jj) - mu/gam*bp2);
        
end

for jj=1:length(Xmesh)

    xvec = 1/2*gvals.*( fphix(Xmesh(jj)-xpos,1 + mu*eta(jj),zpos) );
    zvec = 1/2*gvals.*( fphiz(Xmesh(jj)-xpos,1 + mu*eta(jj),zpos) );
    zveca = 1/2*gvals.*( fphiza(Xmesh(jj)-xpos,1 + mu*eta(jj),zpos) );
    
    phix(jj) = F*sum( xvec );
    phiz(jj) = F*sum( zvec );    
    
    Pv(jj) = phiz(jj) - mu*gam*etax(jj).*phix(jj);
    Ev(jj) = F*sum((xdot.*xvec + gam*zdot.*zveca))-mu*(phix(jj)^2 + phiz(jj)^2)/2;

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

nl2 = (Ev - mu./(1+(gm*etax).^2).*( .5*Q.^2 - mu*gam*(gam*dno+Pv).*term1 +...
                                     phix.*(Q - gam*gm*etax.*dno) - .5*gam*dno.*(2*Pv + gam*dno) + phiz.*(gm*term1+gam*dno) ));

nl1 = fft(dnonl) + fft(Pv)/gam;
nl2 = Dx.*fft(nl2);

% De-alias Bernoulli equation nonlinearity on our way out of this
% subroutine.

nl1(Kc:Kuc) = 0;
nl2(Kc:Kuc) = 0;

nl = [nl1;nl2;xdot;zdot];