function [xpos,zpos,gvals,rvals,Gamma,Nvorts] = initializer_elliptical(Nx,F,av,bv)

xivals = (fliplr(cos(linspace(0,pi/2-1e-1,Nx))))';
dxivals = xivals(2:end) - xivals(1:end-1);
dxivals = [dxivals;dxivals(end)];

%tvals = (linspace(1e-2,2*pi,Nx));

txvals = acos(linspace(-.99999,.99999,Nx/2));
tvals = [txvals (2*pi-txvals)];

c = sqrt(av^2 - bv^2);
xi0 = 1/2*log((av+bv)/(av-bv));

xvals = c*cosh(xi0*xivals)*cos(tvals);
zvals = c*sinh(xi0*xivals)*sin(tvals);
cmat = c*2*pi/Nx*sqrt((cosh(xivals)).^2*ones(1,Nx) - ones(Nx,1)*(cos(tvals)).^2);
csmat = c*sqrt(sqrt((dxivals.*cosh(xivals)).^2*ones(1,Nx) - (dxivals*cos(tvals)).^2));
xpos = xvals(:);
zpos = zvals(:);
rval = .5*(max(max(csmat))+min(min(csmat)));
%rvals = rval*ones(length(xpos),1);
rvals = csmat(:);

Nvorts = length(xpos);
circ = F/Nvorts;

gvals = circ*ones(length(xpos),1);

Gamma = F;
disp('Local Circulation is')
disp(circ)
disp('Number of Starting Vortices is:')
disp(Nvorts)
disp('Froude number is:')
disp(F)
