function [xpos,zpos,gvals,ep,Nvorts] = initializer(Nx,gam,av,omega,zoff)

xvals = linspace(-1,1,Nx);
zvals = 0:2/Nx:1-2/Nx;

disp('Mesh-size is')
disp(2/Nx)
ep = 2*2/Nx;

xpos = [];
zpos = [];
ovals = [];

%cfun = @(x,z) (x/av).^2 + gam^2*((z-zoff)/bv).^2;
%ofun = @(x,z) omega*exp(log(eps)*(cfun(x,z)).^10);
cfun = @(x,z) 1-((x/av).^2 + gam^2*((z-zoff)/av).^2);
ofun = @(x,z) omega*(cfun(x,z)).^3;
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
    ainds = cfun(xvals(jj),zvals) >= 0;
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

rmat = sqrt((xpos*ones(1,Nvorts)-ones(Nvorts,1)*xpos').^2 + (zpos*ones(1,Nvorts)-ones(Nvorts,1)*zpos').^2);
ptmat = 1/ep^2*chi(rmat/ep);
gvals = ptmat\ovals';