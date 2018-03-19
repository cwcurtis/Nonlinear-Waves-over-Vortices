function [xpos,zpos] = init_pos(rows,cols,rows_in,cols_in,n_bdry,Nvorts,xpos1,zpos1,dx,dz,coords,g2,p_shape)

xpos = zeros(Nvorts,1);
zpos = zeros(Nvorts,1);

% Cartesian Grid
if strcmp(coords, 'cartesian')
    for ii = 1:rows
        for jj = 1:cols
            ij = (ii-1)*cols+jj;
            xpos(ij) = xpos1+dx*(jj-(cols+1)/2);
            zpos(ij) = zpos1-dz*(ii-(rows+1)/2);
        end
    end

    %Secondary Grid
    if g2 == 1
        for ii = 1:rows_in
            for jj = 1:cols_in
                ij = rows*cols+(ii-1)*cols_in+jj;
                xpos(ij) = xpos1+dx/d_in*(jj-(cols_in+1)/2);
                zpos(ij) = zpos1-dz/d_in*(ii-(rows_in+1)/2);
            end
        end
    end
    
    % Boundary
    if n_bdry > 0
        i_0 = 1+Nvorts-n_bdry;
        for ii = i_0:Nvorts
            xpos(ii) = xpos1+(cols-1)*dx/2*cos(2*pi*((ii-i_0)/n_bdry));
            zpos(ii) = zpos1+(rows-1)*dz/2*sin(2*pi*((ii-i_0)/n_bdry));
        end
    end

% % % % % % % % % % % % % % % % % %Begin Polar Grid % % % % % % % % % % % % % % % % % % % % %

elseif strcmp(coords, 'polar')
    da = dz/2;
    db = da;        % Change this to get an ellipse
    apos1 = da;
    bpos1 = db;
    thpos1 = pi/2;
    for ii = 1:rows
        for jj = 1:cols
            ij = (ii-1)*cols+jj;
            costh = cos(thpos1+2*pi*(jj-1)/cols);
            sinth = sin(thpos1+2*pi*(jj-1)/cols);
            
            if strcmp(p_shape, 'circle')
                xpos(ij) = xpos1+(apos1+(ii-1)*da)*costh;
                zpos(ij) = zpos1+(bpos1+(ii-1)*db)*sinth;
            end
            
            if strcmp(p_shape, 'cardioid')
                card = 1-cos(2*pi*(jj-1)/cols);
                xpos(ij) = xpos1+(apos1+(ii-1)*da)*costh*card/2;
                zpos(ij) = zpos1+(bpos1+(ii-1)*db)*sinth*card/2;
            end
            
            if strcmp(p_shape, 'rose')
                k = 5;
                rose = cos(k*(2-mod(k,2))*pi*(jj-1)/cols);
                xpos(ij) = xpos1+(apos1+(ii-1)*da)*cos(thpos1+(2-mod(k,2))*pi*(jj-1)/cols)*rose;
                zpos(ij) = zpos1+(bpos1+(ii-1)*db)*sin(thpos1+(2-mod(k,2))*pi*(jj-1)/cols)*rose;
            end
        end
    end
    
    % Boundary
    if n_bdry > 0
        i_0 = 1+Nvorts-n_bdry;
        for ii = i_0:Nvorts
            xpos(ii) = xpos1+cols*da*cos(2*pi*((ii-i_0)/n_bdry));
            zpos(ii) = zpos1+rows*db*sin(2*pi*((ii-i_0)/n_bdry));
        end
    end
end

% % % % % % % % % % % % % % % % % %end polar % % % % % % % % % % % % % % % % % % % % % % % %

%disp([xpos zpos])