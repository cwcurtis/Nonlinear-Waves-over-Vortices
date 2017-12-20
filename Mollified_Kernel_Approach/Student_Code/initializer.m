function [xpos,zpos,gvals,Gamma] = initializer(Nx,omega,cfun,ep)

xvals = linspace(-1,1,Nx);
zvals = linspace(0,1,Nx/2);

circ = omega*(2/Nx)^2;
disp('Mesh-size is')
disp(2/Nx)

disp('Local circulation is')
disp(circ)

xpos = [];
zpos = [];
gvals = [];
Gamma = 0;

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
        vinds = cfun(xvals(jj),zvals(ainds)) <= 1;
        rvec = ones(1,flvds);
        xpos = [xpos xvals(jj)*ones(1,flvds)];
        zpos = [zpos zvals(ainds)];
        gvals = [gvals circ*rvec];
        Gamma = Gamma + circ;        
    end    
end


xpos = xpos.';
zpos = zpos.';
gvals = gvals.';