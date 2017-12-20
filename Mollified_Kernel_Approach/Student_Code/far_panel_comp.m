function Kfar = far_panel_comp(xloc,zloc,xfar,zfar,gfar,gam,mterms)

nparts = length(xloc);
Kfar = zeros(nparts,2);
pg = pi*gam;
pvls = zeros(nparts,mterms);
jvals = 1:mterms;

xeval = xloc*jvals;
xevalf = xfar*jvals;

Zmat = zloc*ones(1,length(zfar)) - ones(nparts,1)*zfar';
Zplus = Zmat > 0;
Zminus = Zmat < 0;
Zequal = Zmat == 0;


exeval = exp(1i*pi*xeval);
exevalf = exp(-1i*pi*xevalf);

zplus = zeros(nparts,mterms);
znegs = zeros(nparts,mterms);
zequs = zeros(nparts,mterms);

for ll=1:nparts
   
    if ll==1        
        for jj = 1:mterms               
            zpnew = sum(gfar(Zplus(ll,:)).*exevalf(Zplus(ll,:),jj).*exp(pg*jj*zfar(Zplus(ll,:))));
            znnew = sum(gfar(Zminus(ll,:)).*exevalf(Zminus(ll,:),jj).*exp(-pg*jj*zfar(Zminus(ll,:))));
            zenew = sum(gfar(Zequal(ll,:)).*exevalf(Zequal(ll,:),jj));

            pvls(ll,jj) = zpnew*exp(-pg*zloc(ll)) + znnew*exp(pg*zloc(ll)) + zenew;        
            
            zplus(1,jj) = zpnew;
            znegs(1,jj) = znnew;
            zequs(1,jj) = zenew;            
        end
    else
        flagp = 0;
        jjp = 1;
        while jjp < ll          
            if isequal(Zplus(jjp,:),Zplus(ll,:))
                flagp = 1;
            else
                jjp = jjp+1;
            end           
        end
        
        flagn = 0;
        jjn = 1;
        while jjn < ll           
            if isequal(Zminus(jjn,:),Zminus(ll,:))
                flagn = 1;
            else
                jjn = jjn+1;
            end           
        end
        
        flage = 0;
        jje = 1;
        while jje < ll           
            if isequal(Zequal(jje,:),Zequal(ll,:))
                flage = 1;
            else
                jje = jje+1;
            end           
        end
        
        for jj = 1:mterms               
            if flagp == 0
                zpnew = sum(gfar(Zplus(ll,:)).*exevalf(Zplus(ll,:),jj).*exp(pg*jj*zfar(Zplus(ll,:))));
                zplus(ll,jj) = zpnew;
            else
                zpnew = zplus(jjp,jj);
            end
            
            if flagn == 0
                znnew = sum(gfar(Zminus(ll,:)).*exevalf(Zminus(ll,:),jj).*exp(-pg*jj*zfar(Zminus(ll,:))));
                znegs(ll,jj) = znnew;
            else
                znnew = znegs(jjn,jj);
            end
            
            if flage == 0            
                zenew = sum(gfar(Zequal(ll,:)).*exevalf(Zequal(ll,:),jj));
                zequs(ll,jj) = zenew;
            else
                zenew = zequs(jje,jj);
            end

            pvls(ll,jj) = zpnew*exp(-pg*zloc(ll)) + znnew*exp(pg*zloc(ll)) + zenew;                
        end       
    end   
end

Kvals = sum(exeval.*pvls,2);
Kfar(:,1) = 2*real(Kvals);
Kfar(:,2) = -2*imag(Kvals);

