function gvals = set_gvals(rows,cols,rows_in,cols_in,Nvorts,c1,c2,c3,c4,c5,pattern,g2)

gvals = zeros(Nvorts,1);

% Automatically assign vortex strengths
for ii = 1:rows
    for jj = 1:cols
        ij = (ii-1)*cols+jj;
        
        if strcmp(pattern, 'uniform')
            gvals(ij) = c1;
        end
        
%        if mod(jj,2) == 0
%            gvals(ij) = -gvals(ij);
%        end
        
        a = (c1+c2)/2;
        b = (c1-c2)/2;
        
        if strcmp(pattern, 'hor_stripes')
            gvals(ij) = a+b*(-1)^ii;
        end
        
        if strcmp(pattern, 'ver_stripes')
            gvals(ij) = a+b*(-1)^jj;
        end
        
        if strcmp(pattern, 'checkerboard')
            gvals(ij) = a+b*(-1)^(ii+jj);
        end
        
        if strcmp(pattern, 'clustered')
            gvals(ij) = (jj < cols/3)*c1+(jj > cols*2/3)*c2;
        end
        
        if strcmp(pattern, 'gradient')
            a1 = (rows-ii)*(cols-jj)*c1;
            a2 = (rows-ii)*(jj-1)*c2;
            a3 = (ii-1)*(cols-jj)*c3;
            a4 = (ii-1)*(jj-1)*c4;
            a5 = (rows-1)*(cols-1);
            gvals(ij) = (a1+a2+a3+a4)/a5;
        end
    end
end

% Secondary Grid gvals
if g2 == 1
    for ii = 1+rows*cols:rows_in*cols_in+rows*cols
        gvals(ii) = c5;
    end
end