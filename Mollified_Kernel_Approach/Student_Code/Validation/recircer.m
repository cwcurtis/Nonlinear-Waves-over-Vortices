function [xpud,zpud,gvud] = recircer(gvals,xpos,zpos,Nx)

dx = 2/Nx;

xl = min(xpos);
xr = max(xpos);
zb = min(zpos);
zt = max(zpos);

xgl = -1 + floor((xl+1)/dx)*dx;
zgb = -.5 + floor((zb+.5)/dx)*dx;

hcnt = ceil((xr+1)/dx) - floor((xl+1)/dx);
vcnt = ceil(zt/dx) - floor(zb/dx);

xpud = [];
zpud = [];
gvud = [];
disp('Circulation going in')
disp(sum(gvals))

for jj = 1:hcnt + 1
    for kk = 1:vcnt + 1
        hpc = xgl + (jj-1)*dx;
        vpc = zgb + (kk-1)*dx;
        
        hindsl = xpos >= hpc - 2*dx;
        hindsr = xpos <= hpc + 2*dx;
        vindsb = zpos >= vpc - 2*dx;
        vindst = zpos <= vpc + 2*dx;
        linds = logical(hindsl.*hindsr.*vindsb.*vindst);
        
        if length(zpos(linds))>1
           xpud = [xpud hpc];
           zpud = [zpud vpc];
           glocs = gvals(linds);
           
           zdif = abs(zpos(linds) - vpc)/dx;
           xdif = abs(xpos(linds) - hpc)/dx;
           zdif2 = zdif.*zdif;
           zdif3 = zdif2.*zdif;
           xdif2 = xdif.*xdif;
           xdif3 = xdif2.*xdif;
                      
           zlto = zdif<=1;
           zgto = zdif>1;
           xlto = xdif<=1;
           xgto = xdif>1;
           
           Rg1 = logical(zlto.*xlto);
           Rg2 = logical(zlto.*xgto);
           Rg3 = logical(zgto.*xlto);
           Rg4 = logical(zgto.*xgto);
           
           if(length(Rg1)>=1)
               gRg1 = glocs(Rg1).*(1 - 5/2*xdif2(Rg1)+3/2*xdif3(Rg1)).*(1 - 5/2*zdif2(Rg1)+3/2*zdif3(Rg1));
           else
               gRg1 = 0;
           end
           
           if(length(Rg2)>=1)
               gRg2 = glocs(Rg2).*(2-4*xdif(Rg2)+5/2*xdif2(Rg2)-xdif3(Rg2)/2).*(1 - 5/2*zdif2(Rg2)+3/2*zdif3(Rg2));
           else
               gRg2 = 0;
           end
           
           if(length(Rg3)>=1)
               gRg3 = glocs(Rg3).*(1 - 5/2*xdif2(Rg3)+3/2*xdif3(Rg3)).*(2-4*zdif(Rg3)+5/2*zdif2(Rg3)-zdif3(Rg3)/2);
           else
               gRg3 = 0;
           end
           
           if(length(Rg4)>=1)
               gRg4 = glocs(Rg4).*(2-4*xdif(Rg4)+5/2*xdif2(Rg4)-xdif3(Rg4)/2).*(2-4*zdif(Rg4)+5/2*zdif2(Rg4)-zdif3(Rg4)/2);
           else
               gRg4 = 0;
           end
           
           gup = sum(gRg1) + sum(gRg2) + sum(gRg3) + sum(gRg4);
           gvud = [gvud gup];
        end
    end
end

disp('Circulation coming out')
disp(sum(gvud))
%{
clf
hold on
plot(1:length(gvals),log10(abs(gvals)))
plot(1:length(gvud),log10(abs(gvud)))
hold off
pause
%}
xpud = xpud';
zpud = zpud';
gvud = gvud';

