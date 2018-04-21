function [xpud,zpud,gvud] = recircer_bndry(gvals,xpos,zpos,Nx)

dx = 2/Nx;

xmesh = -1:dx:1-dx;
zmesh = dx:dx:1-dx;
zul = ceil(.9/dx)*dx;
zlbnds = zpos < 1.5*dx;
zubnds = zpos > zul-1.5*dx;
zlow = zpos(zlbnds);
zhigh = zpos(zubnds);
zlbndsdx = zpos < 2*dx;
zubndsdx = zpos > zul-2*dx;
zlowdx = zpos(zlbndsdx);
zhighdx = zpos(zubndsdx);

xl = min(xpos);
xr = max(xpos);
zrem = zpos(logical(1-zlbnds.*zubnds));
zb = min(zrem);
zt = max(zrem);

zlflag = 1;
if isempty(zlow)
    zlflag = 0;
    lsect = 0;
end

zuflag = 1;
if isempty(zhigh)
    zuflag = 0;
    usect = 0;
end

jjl = max(floor((xl+1)/dx)-1,1);
jju = min(ceil((xr+1)/dx)+1,Nx);

kkl = max(floor(zb/dx)-1,2);
kku = min(ceil(zt/dx)+1,ceil(.9/dx)-1);

xpud = [];
zpud = [];
gvud = [];
for jj = jjl:jju
    hpc = xmesh(jj);
    hinds = abs(xpos-hpc) <= 3/2*dx;
    if zlflag == 1
        llinds = logical(hinds.*zlbnds);
        glocsl = gvals(llinds);           
        lsect = sum(llinds);
        xdifl = (xpos(llinds) - hpc)/dx;
        xdifl2 = xdifl.*xdifl;
        xltol = abs(xdifl)<=1/2;
        xgtol = logical(1-xltol);
        glsll = glocsl(xltol);
        %xdsll = xdifl(xltol);
        xdssqll = xdifl2(xltol);                                    
        glsgl = glocsl(xgtol);
        xdsgl = xdifl(xgtol);
        xdssqgl = xdifl2(xgtol);                                    
    end
    if zuflag == 1
        ulinds = logical(hinds.*zubnds);       
        glocsu = gvals(ulinds);           
        usect = sum(ulinds);
        xdifu = (xpos(ulinds) - hpc)/dx;
        xdifu2 = xdifu.*xdifu;
        xltou = abs(xdifu)<=1/2;
        xgtou = logical(1-xltou);
        glslu = glocsu(xltou);
        %xdslu = xdifu(xltou);
        xdssqlu = xdifu2(xltou);                                    
        glsgu = glocsu(xgtou);
        xdsgu = xdifu(xgtou);
        xdssqgu = xdifu2(xgtou);                                    
    end
    for kk = kkl:kku
        vpc = zmesh(kk);
        vinds = abs(zrem-vpc) <= 3/2*dx;
        linds = logical(hinds.*vinds);
        
        if ~isempty(zrem(linds))
           xpud = [xpud hpc];
           zpud = [zpud vpc];
           glocs = gvals(linds);
           
           zdif = abs(zrem(linds) - vpc)/dx;
           xdif = abs(xpos(linds) - hpc)/dx;
           zdif2 = zdif.*zdif;
           xdif2 = xdif.*xdif;
                      
           zlto = zdif<=1/2;
           zgto = logical(1-zlto);
           xlto = xdif<=1/2;
           xgto = logical(1-xlto);
           
           Rg1 = logical(zlto.*xlto);
           Rg2 = logical(zlto.*xgto);
           Rg3 = logical(zgto.*xlto);
           Rg4 = logical(zgto.*xgto);
           
           if length(Rg1)>1
               gRg1 = glocs(Rg1).*(1 - xdif2(Rg1)).*(1 - zdif2(Rg1));
           else
               gRg1 = 0;
           end
           
           if length(Rg2)>1
               gRg2 = glocs(Rg2).*(1 - 3*xdif(Rg2)/2 + xdif2(Rg2)/2).*(1 - zdif2(Rg2));
           else
               gRg2 = 0;
           end
           
           if length(Rg3)>1
               gRg3 = glocs(Rg3).*(1 - xdif2(Rg3)).*(1 - 3*zdif(Rg3)/2 + zdif2(Rg3)/2);
           else
               gRg3 = 0;
           end
           
           if length(Rg4)>1
               gRg4 = glocs(Rg4).*(1 - 3*xdif(Rg4)/2 + xdif2(Rg4)/2 ).*(1 - 3*zdif(Rg4)/2 + zdif2(Rg4)/2 );
           else
               gRg4 = 0;
           end  
           
           gup = sum(gRg1) + sum(gRg2) + sum(gRg3) + sum(gRg4);
           
           if lsect > 0
                if vpc == 2*dx 
                    zdif2dx = abs(zlow - 2*dx)/dx;
                    zdif2dx2 = zdif2dx.*zdif2dx;                                
                    if length(xltol)>1
                        gRg12 = glsll.*(1 - xdssqll).*(1 - zdif2dx2(xltol));                       
                    else
                        gRg12 = 0;                
                    end           
                    if length(xgtol)>1
                        gRg22 = glsgl.*(1 - 3*xdsgl/2 + xdssqgl/2).*(1 - zdif2dx2(xgtol));               
                    else
                        gRg22 = 0;                
                    end                     
                    gup = gup + sum(gRg12) + sum(gRg22);
                end
           
                if vpc == 3*dx
                    zdif3dx = abs(zlow - 3*dx)/dx;
                    zdif3dx2 = zdif3dx.*zdif3dx;                        
                    if length(xltol)>1
                        gRg13 = glsll.*(1 - xdssqll).*(1 - 3/2*zdif3dx(xltol) + 1/2*zdif3dx2(xltol));                       
                    else
                        gRg13 = 0;                
                    end           
                    if length(xgtol)>1
                        gRg23 = glsgl.*(1 - 3*xdsgl/2 + xdssqgl/2).*(1 - 3/2*zdif3dx(xgtol) + 1/2*zdif3dx2(xgtol));               
                    else
                        gRg23 = 0;                
                    end     
                    gup = gup + sum(gRg13) + sum(gRg23);
                end
           end           
           
           if usect > 0
                if vpc == zul-2*dx
                    zdif2dx = abs(zhigh - (zul-2*dx))/dx;
                    zdif2dx2 = zdif2dx.*zdif2dx;           
                    if length(xltou)>1
                        gRg12 = glslu.*(1 - xdssqlu).*(1 - zdif2dx2(xltou));                       
                    else
                        gRg12 = 0;                
                    end          
                    if length(xgtou)>1
                        gRg22 = glsgu.*(1 - 3*xdsgu/2 + xdssqgu/2).*(1 - zdif2dx2(xgtou));               
                    else
                        gRg22 = 0;                
                    end                           
                    gup = gup + sum(gRg12) + sum(gRg22);
                end
           
                if vpc == zul-3*dx
                    zdif3dx = abs(zhigh - (zul-3*dx))/dx;
                    zdif3dx2 = zdif3dx.*zdif3dx;       
                    if length(xltou)>1
                        gRg13 = glslu.*(1 - xdssqlu).*(1 - 3/2*zdif3dx(xltou) + 1/2*zdif3dx2(xltou));                       
                    else
                        gRg13 = 0;                
                    end           
                    if length(xgtou)>1
                        gRg23 = glsgu.*(1 - 3*xdsgu/2 + xdssqgu/2).*(1 - 3/2*zdif3dx(xgtou) + 1/2*zdif3dx2(xgtou));               
                    else
                        gRg23 = 0;                
                    end     
                    gup = gup + sum(gRg13) + sum(gRg23);
                end
           end           
           
           gvud = [gvud gup];
        end
    end
    
    if lsect>0
       xpud = [xpud hpc];
       zpud = [zpud dx];            
       zdifdx = abs(zlowdx - dx)/dx;
       zdifdx2 = zdifdx.*zdifdx;            
       if length(xltol)>1
          gRg11 = glsll.*(1 - xdssqll).*(1 - 3/2*zdifdx(xltol) + 1/2*zdifdx2(xltol));                       
       else
          gRg11 = 0;                
       end           
       if length(xgtol)>1
          gRg21 = glsgl.*(1 - 3*xdsgl/2 + xdssqgl/2).*(1 - 3/2*zdifdx(xgtol) + 1/2*zdifdx2(xgtol));               
       else
          gRg21 = 0;                
       end            
       gup1 = sum(gRg11) + sum(gRg21);
       gvud = [gvud gup1];            
    end
    
    if usect>0
       xpud = [xpud hpc];
       zpud = [zpud zul-dx];           
       zdifdx = abs(zhighdx - (zul-dx))/dx;
       zdifdx2 = zdifdx.*zdifdx;            
       if length(xltou)>1
          gRg11 = glslu.*(1 - xdssqlu).*(1 - 3/2*zdifdx(xltou) + 1/2*zdifdx2(xltou));                
       else
          gRg11 = 0;                
       end           
       if length(xgtou)>1
          gRg21 = glsgu.*(1 - 3*xdsgu/2 + xdssqgu/2).*(1 - 3/2*zdifdx(xgtou) + 1/2*zdifdx2(xgtou));                
       else
          gRg21 = 0;                
       end           
       gup1 = sum(gRg11) + sum(gRg21);            
       gvud = [gvud gup1];                  
    end
    
end
xpud = xpud';
zpud = zpud';
gvud = gvud';