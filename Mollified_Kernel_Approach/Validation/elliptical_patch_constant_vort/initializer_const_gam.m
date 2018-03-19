function [xpos,zpos,gvals,ep,Gamma,Nvorts] = initializer_const_gam(Nx,F,av,bv)

xvals = linspace(-1,1,Nx);
zvals = -1/2:2/Nx:(1/2-2/Nx);
dx = 2/Nx;
disp('Mesh-size is')
disp(dx)
ep = dx;

xpos = [];
zpos = [];
gvals = [];

cfun = @(x,z) (x/av).^2 + (z/bv).^2;

%% In the off chance I want ghost particles later. 
%{
for jj=1:Nx    
    ainds = cfun(xvals(jj),zvals) <= 1+10*ep;
    flvds = length(zvals(ainds));
    if flvds>=1
        vinds = cfun(xvals(jj),zvals(ainds)) <= 1;
        rvec = zeros(1,flvds);
        if length(vinds)>=1
            rvec(vinds) = 1;
        end
        xpos = [xpos xvals(jj)*ones(1,flvds)];
        zpos = [zpos zvals(ainds)];
        gvals = [gvals rvec];        
    end    
end
%}


for jj=1:Nx    
    ainds = cfun(xvals(jj),zvals) <= 1;
    flvds = length(zvals(ainds));
    if flvds>=1
        xpos = [xpos xvals(jj)*ones(1,flvds)];
        zpos = [zpos zvals(ainds)];
        gvals = [gvals ones(1,flvds)];
    end    
end


Nvorts = length(gvals);
circ = F/Nvorts;
Gamma = F;
disp('Number of Starting Vortices is:')
disp(Nvorts)
disp('Kernel Width is:')
disp(ep)
xpos = xpos.';
zpos = zpos.';
gvals = circ*gvals.';
