function mp = mintper(xi,dx,ntrunc)

phi = @(k) (9 - 12*cos(k) + k.*sin(2*k)-2*k.*sin(k) + 3*cos(2*k))/k.^4;

phim = repmat(pi*dx*phi(1:ntrunc),length(dx),1);

mp = dx*(1/2 + 2*sum(cos(pi*xi*(1:ntrunc)).*phim,2));