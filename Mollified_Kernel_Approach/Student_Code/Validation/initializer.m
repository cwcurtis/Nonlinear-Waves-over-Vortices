function [xpos,zpos,gvals,Gamma] = initializer(Nx,omega,cfun,gam,Rv,zoff)

xvals = linspace(-1,1,Nx);
zvals = linspace(-1,2,3*Nx/2);

disp('Mesh-size is')
disp(2/Nx)

%disp('Local circulation is')
%disp(circ)

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
    ainds = cfun(xvals(jj),zvals) <= 1;
    flvds = length(zvals(ainds));
    if flvds>=1
        vinds = cfun(xvals(jj),zvals(ainds)) <= 1;
        rvec = ones(1,flvds);
        xpos = [xpos xvals(jj)*ones(1,flvds)];
        zpos = [zpos zvals(ainds)];
        circ = omega*(Rv^2 - (xvals(jj))^2 - gam^2*(zvals(ainds)-zoff).^2).^3*(2/Nx)^2;
        gvals = [gvals circ];          
    end    
end

Gamma = pi*omega*Rv^8/(4*gam);
disp('Total circulation is')
disp(Gamma)
xpos = xpos.';
zpos = zpos.';
gvals = gvals.';