function [xpos,zpos,gvals,Nvorts,ep] = initializer(Nx,om,F,gam,cfun,Rv)

xvals = linspace(-1,1,Nx);
zvals = -1/2:2/Nx:(1/2-2/Nx);
dx = 2/Nx;
disp('Mesh-size is')
disp(dx)

xpos = [];
zpos = [];
gvals = [];

omfunc = @(xj,yj) om/(F*gam)*dx^2*(Rv^2 - xj.^2 - yj.^2);

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
    ainds = cfun(xvals(jj),zvals,Rv) <= 1;
    flvds = length(zvals(ainds));
    if flvds>=1
        xpos = [xpos xvals(jj)*ones(1,flvds)];
        zpos = [zpos zvals(ainds)];
        gvals = [gvals omfunc(xvals(jj),zvals(ainds))];
    end    
end

Nvorts = length(gvals);
ep = 2/Nx;

disp('Number of Starting Vortices is:')
disp(Nvorts)
disp('Kernel Width is:')
disp(ep)
xpos = xpos.';
zpos = zpos.';
gvals = gvals.';