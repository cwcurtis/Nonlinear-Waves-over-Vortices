function [xpos,zpos,gvals,ep,Nvorts] = initializer(Nx,omega,gam,Rv)

xvals = linspace(-1,1,Nx);
zvals = -1/2:2/Nx:(1/2-2/Nx);
dx = 2/Nx;
disp('Mesh-size is')
disp(dx)
ep = .95*dx;

xpos = [];
zpos = [];
ovals = [];

ofun = @(xj,yj) omega*(1. - (xj.^2 + yj.^2)/Rv.^2).^3;
cfun = @(x,z) (x.^2 + z.^2)/Rv^2;
chi = @(r) 1./pi.*(2.*exp(-(r.^2))-.5*exp(-(r/sqrt(2)).^2));
    
%% In the off chance I want ghost particles later. 
%{
for jj=1:Nx    
    ainds = cfun(xvals(jj),zvals) <= 1+2*ep;
    flvds = length(zvals(ainds));
    if flvds>=1
        vinds = cfun(xvals(jj),zvals(ainds)) <= 1;
        rvec = zeros(1,flvds);
        if length(vinds)>=1
            rvec(vinds) = circ;
        end
        xpos = [xpos xvals(jj)*ones(1,flvds)];
        zpos = [zpos zvals(ainds)];
        gvals = [gvals rvec];
        Gamma = Gamma + circ;        
    end    
end
%}

for jj=1:Nx    
    ainds = cfun(xvals(jj),zvals) <= 1;
    flvds = length(zvals(ainds));
    if flvds>=1
        xpos = [xpos xvals(jj)*ones(1,flvds)];
        zpos = [zpos zvals(ainds)];
        ovals = [ovals ofun(xvals(jj),zvals(ainds))/gam];
    end    
end

Nvorts = length(ovals);
disp('Number of Starting Vortices is:')
disp(Nvorts)
disp('Kernel Width is:')
disp(ep)
xpos = xpos.';
zpos = zpos.';

rmat = sqrt((repmat(xpos,1,Nvorts)-repmat(xpos',Nvorts,1)).^2 + (repmat(zpos,1,Nvorts)-repmat(zpos',Nvorts,1)).^2);
ptmat = 1/ep^2*chi(rmat/ep);
gvals = ptmat\ovals';