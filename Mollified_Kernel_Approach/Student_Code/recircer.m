function [xpud,zpud,gvud] = recircer(gvals,xpos,zpos,Nx)

dx = 2/Nx;

xl = min(xpos);
xr = max(xpos);
zb = min(zpos);
zt = max(zpos);

xgl = floor(xl/dx)*dx;
zgb = floor(zb/dx)*dx;

hcnt = ceil(xr/dx) - floor(xl/dx);
vcnt = ceil(zt/dx) - floor(zb/dx);

xpud = [];
zpud = [];
gvud = [];

for jj = 1:hcnt + 1
    for kk = 1:vcnt + 1
        hpc = xgl + (jj-1)*dx;
        vpc = zgb + (kk-1)*dx;
        
        hindsl = xpos >= hpc - 2*dx;
        hindsr = xpos <= hpc + 2*dx;
        vindsb = zpos >= vpc - 2*dx;
        vindst = zpos <= vpc + 2*dx;
        linds = logical(hindsl.*hindsr.*vindsb.*vindst);
        
        if length(zpos(linds))>=1
           xpud = [xpud hpc];
           zpud = [zpud vpc];
           glocs = gvals(linds);
           
           zdif = abs(zpos(linds) - vpc)/dx;
           xdif = abs(xpos(linds) - hpc)/dx;
           zdif2 = zdif.*zdif;
           zdif3 = zdif2.*zdif;
                      
           zlto = zdif<=1;
           zgto = zdif>1;
           zvlto = zdif(zlto);
           zvgto = zdif(zgto);
           
           xint = mintper(xdif,dx,10);
           
           if(length(zvlto)>=1)
               glto = glocs(zlto).*xint(zlto).*(1 - 5/2*zdif2(zlto)+3/2*zdif3(zlto));
           else
               glto = 0;
           end
           
           if(length(zvgto)>=1)
               ggto = glocs(zgto).*xint(zgto).*(2-4*zdif(zgto)+5/2*zdif2(zgto)-zdif3(zgto)/2);
           else
               ggto = 0;
           end
           
           gup = sum(glto) + sum(ggto);
           gvud = [gvud gup];
        end
    end
end

