function interp_test(Nx,Nxp)

dx = 2/Nx;

xpos = linspace(-.95,.95,Nxp)';

disp(1/Nxp)

xl = min(xpos);
xr = max(xpos);

xgl = floor(xl/dx)*dx;
hcnt = ceil(xr/dx) - floor(xl/dx);

for jj = 1:hcnt + 1
        hpc = xgl + (jj-1)*dx;
        
        hindsl = xpos >= hpc - 2*dx;
        hindsr = xpos <= hpc + 2*dx;
        linds = logical(hindsl.*hindsr);
        
        if length(xpos(linds))>=1
           
           xdif = abs(xpos(linds)-hpc)/dx;
           xint = mintper(xdif,dx,400);
          
           plot(xpos(linds),xint)
           pause
           
        end    
end

