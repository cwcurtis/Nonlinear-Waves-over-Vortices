function [xpos,zpos,gvals,Gamma,Nvorts] = initializer(Nx,F,cfun,gam,av,bv,zoff)

xvals = linspace(-1,1,Nx);
zvals = 0:2/Nx:1-2/Nx;

disp('Mesh-size is')
disp(2/Nx)

xpos = [];
zpos = [];
gvals = [];

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
    ainds = cfun(xvals(jj),zvals,gam,av,bv,zoff) <= 1;
    flvds = length(zvals(ainds));
    if flvds>=1
        vinds = cfun(xvals(jj),zvals(ainds),gam,av,bv,zoff) <= 1;
        rvec = ones(1,flvds);
        xpos = [xpos xvals(jj)*ones(1,flvds)];
        zpos = [zpos zvals(ainds)];
        gvals = [gvals rvec];
    end    
end

Nvorts = length(gvals);
circ = F/Nvorts;
Gamma = F;
disp('Local Circulation is')
disp(circ)
disp('Number of Starting Vortices is:')
disp(Nvorts)
disp('Froude number is')
disp(F)
        
xpos = xpos.';
zpos = zpos.';
gvals = circ*gvals.';