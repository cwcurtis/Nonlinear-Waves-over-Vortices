function [xpud,zpud,gvud] = recircer_bndry(gvals,xpos,zpos,Nx)

dx = 2/Nx;

xl = min(xpos);
xr = max(xpos);
zb = min(zpos);
zt = max(zpos);

xgl = -1  + floor((xl+1)/dx)*dx;
zgb = floor(zb/dx)*dx;

hcnt = ceil((xr+1)/dx) - floor((xl+1)/dx);
vcnt = ceil(zt/dx) - floor(zb/dx);

xpud = [];
zpud = [];
gvud = [];

for jj = 1:hcnt + 1
    for kk = 1:vcnt + 1
        hpc = xgl + (jj-1)*dx;
        vpc = zgb + (kk-1)*dx;
        hindsl = xpos >= hpc - 3/2*dx;
        hindsr = xpos <= hpc + 3/2*dx;
        vindsb = zpos >= vpc - 3/2*dx;
        vindst = zpos <= vpc + 3/2*dx;
        linds = logical(hindsl.*hindsr.*vindsb.*vindst);
        
        if length(zpos(linds))> 1
           xpud = [xpud hpc];
           zpud = [zpud vpc];
           glocs = gvals(linds);
           
           zdif = abs(zpos(linds) - vpc)/dx;
           xdif = abs(xpos(linds) - hpc)/dx;
           zdif2 = zdif.*zdif;
           xdif2 = xdif.*xdif;
                      
           zlto = zdif<=1/2;
           zgto = zdif>1/2;
           xlto = xdif<=1/2;
           xgto = xdif>1/2;
           
           Rg1 = logical(zlto.*xlto);
           Rg2 = logical(zlto.*xgto);
           Rg3 = logical(zgto.*xlto);
           Rg4 = logical(zgto.*xgto);
           
           if(length(Rg1)>1)
               gRg1 = glocs(Rg1).*(1 - xdif2(Rg1)).*(1 - zdif2(Rg1));
           else
               gRg1 = 0;
           end
           
           if(length(Rg2)>1)
               gRg2 = glocs(Rg2).*(1 - 3*xdif(Rg2)/2 + xdif2(Rg2)/2).*(1 - zdif2(Rg2));
           else
               gRg2 = 0;
           end
           
           if(length(Rg3)>1)
               gRg3 = glocs(Rg3).*(1 - xdif2(Rg3)).*(1 - 3*zdif(Rg3)/2 + zdif2(Rg3)/2);
           else
               gRg3 = 0;
           end
           
           if(length(Rg4)>1)
               gRg4 = glocs(Rg4).*(1 - 3*xdif(Rg4)/2 + xdif2(Rg4)/2 ).*(1 - 3*zdif(Rg4)/2 + zdif2(Rg4)/2 );
           else
               gRg4 = 0;
           end
           
           gup = sum(gRg1) + sum(gRg2) + sum(gRg3) + sum(gRg4);
           gvud = [gvud gup];
        end
    end
end

xpud = xpud';
zpud = zpud';
gvud = gvud';

