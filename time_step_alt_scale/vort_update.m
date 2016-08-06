function rhs = vort_update(Xmesh,gam,mu,eta,Q,xpos,zpos,gvals,dno)

KT = length(Xmesh);
xdot = zeros(length(xpos),1);
zdot = zeros(length(xpos),1);

fpsix = @(x,z,zj) -.25*sin(pi*x).*( 1./( cosh(pi*gam*(z-zj)) - cos(pi*x) ) + 1./( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );
fpsiz = @(x,z,zj) -.25*( sinh(pi*gam*(z-zj))./( cosh(pi*gam*(z-zj)) - cos(pi*x) ) + sinh(pi*gam*(z+zj))./( cosh(pi*gam*(z+zj)) - cos(pi*x) ) );

for jj=1:length(xpos)
   
    for ll=1:length(xpos)
       
        if(ll ~= jj)
            fac = ( cosh(pi*gam*(zpos(jj)-zpos(ll))) - cos(pi*(xpos(jj)-xpos(ll))) ).*( cosh(pi*gam*(zpos(jj)+zpos(ll))) - cos(pi*(xpos(jj)-xpos(ll))) );
            xdot(jj) = xdot(jj) + gvals(ll)*sinh(pi*gam*zpos(ll))*(cosh(pi*gam*zpos(ll))-cosh(pi*gam*zpos(jj))*cos(xpos(jj)-xpos(ll)))/fac;
            zdot(jj) = zdot(jj) + gvals(ll)*sin(xpos(jj)-xpos(ll))*sinh(pi*gam*zpos(ll))/fac;
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

rhs = [xdot;zdot];