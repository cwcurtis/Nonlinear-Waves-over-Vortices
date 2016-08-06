function [Pv,nl2] = force_terms(Xmesh,gam,mu,eta,etax,Q,xpos,zpos,gvals,dno)

KT = length(Xmesh);
K = KT/2;
Kc = floor(2*K/3);
Kuc = 2*K-Kc+1;
Kc = Kc + 1;
Dx = pi*1i*[0:K -K+1:-1]'; 

xdot = zeros(length(xpos),1);
zdot = zeros(length(xpos),1);

Pv = zeros(length(Xmesh),1);
Ev = zeros(length(Xmesh),1);

phix = zeros(length(Xmesh),1);
phiz = zeros(length(Xmesh),1);

fphiz = @(x,z,zj) sin(pi*x).*sinh(pi*gam*z).*sinh(pi*gam*zj)./( ( cosh(pi*gam*(z-zj)) - cos(pi*x) ).*( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );
fphix = @(x,z,zj) sinh(pi*gam*zj).*(cosh(pi*gam*zj) - cosh(pi*gam*z).*cos(pi*x))./( (cosh(pi*gam*(z-zj)) - cos(pi*x)).*( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );

fpsix = @(x,z,zj) -.25*sin(pi*x).*( 1./( cosh(pi*gam*(z-zj)) - cos(pi*x) ) + 1./( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );
fpsiz = @(x,z,zj) -.25*( sinh(pi*gam*(z-zj))./( cosh(pi*gam*(z-zj)) - cos(pi*x) ) + sinh(pi*gam*(z+zj))./( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );

for jj=1:length(xpos)
   
    for ll=1:length(xpos)
       
        if(ll ~= jj)
            fac = ( cosh(pi*gam*(zpos(jj)-zpos(ll))) - cos(pi*(xpos(jj)-xpos(ll))) ).*( cosh(pi*gam*(zpos(jj)+zpos(ll))) - cos(pi*(xpos(jj)-xpos(ll))) );
            xdot(jj) = xdot(jj) + gvals(ll)*sinh(pi*gam*zpos(ll))*(cosh(pi*gam*zpos(ll))-cosh(pi*gam*zpos(jj))*cos(pi*(xpos(jj)-xpos(ll))))/fac;
            zdot(jj) = zdot(jj) + gvals(ll)*sin(pi*(xpos(jj)-xpos(ll)))*sinh(pi*gam*zpos(ll))/fac;
        end
        
    end
    
    xdot(jj) = mu/4*( 2*xdot(jj) + gvals(jj)/(tanh(gam*pi*zpos(jj))) );
    zdot(jj) = mu/(2*gam)*sinh(gam*pi*zpos(jj))*zdot(jj);    
    
    psix = fpsix(Xmesh-xpos(jj),1+mu*eta,zpos(jj));
    psiz = fpsiz(Xmesh-xpos(jj),1+mu*eta,zpos(jj));
    
    bp1 = 2*sum(-gam*dno.*psix - Q.*psiz)/KT;
    bp2 = 2*sum(-gam*dno.*psiz + Q.*psix)/KT; 
    
    xdot(jj) = xdot(jj) + mu*bp1;
    zdot(jj) = zdot(jj) + mu/gam*bp2;
    
end

for jj=1:length(Xmesh)

    xvec = 1/2*gvals.*( fphix(Xmesh(jj)-xpos,1 + mu*eta(jj),zpos) );
    zvec = 1/2*gvals.*( fphiz(Xmesh(jj)-xpos,1 + mu*eta(jj),zpos) );
    
    phix(jj) = sum( xvec );
    phiz(jj) = sum( zvec );    
    
    Pv(jj) = phiz(jj) - mu*gam*etax(jj).*phix(jj);
    Ev(jj) = sum(gvals.*(xdot.*xvec + gam*zdot.*zvec)) - mu*(phix(jj)^2 + phiz(jj)^2)/2;
    
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

etat = dno + Pv/gam;

nl2 = (Ev - mu./(1+(gm*etax).^2).*( 1/2*Q.^2 - mu*gam^2*etat.*term1 +...
                                     phix.*(Q + gm*etax.*(Pv-gam*etat)) + .5*(Pv.^2-gam^2*etat.^2) + phiz.*(gm*term1-(Pv - gam*etat)) ));

nl2 = Dx.*fft(nl2);

Pv = fft(Pv);

% De-alias Bernoulli equation nonlinearity on our way out of this
% subroutine.

nl2(Kc:Kuc) = 0;