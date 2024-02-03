function ukp1 = dmm(uk, dt, rval, omega, gvals, Nvorts)

ukp10 = uk + dt*force_terms_on_molly_non_periodic(rval,omega,uk,gvals,Nvorts);

fopt = @(u) coords(u, uk, dt, gvals, Nvorts);

options = optimoptions('fsolve','Display','none');

ukp1 = fsolve(fopt, ukp10, options);

end

function system = coords(ukp1, uk, dt, gvals, Nvorts)

xk = uk(1:Nvorts);
yk = uk(Nvorts+1:2*Nvorts);

xkp1 = ukp1(1:Nvorts);
ykp1 = ukp1(Nvorts+1:2*Nvorts);

i2p = 1./(2.*pi);

dxk = bsxfun(@minus,xk,xk');
dyk = bsxfun(@minus,yk,yk');

dxkp1 = bsxfun(@minus,xkp1,xkp1');
dykp1 = bsxfun(@minus,ykp1,ykp1');
    
xmid = .5*(dxk + dxkp1);
ymid = .5*(dyk + dykp1);

rk = dxk.^2 + dyk.^2;
rkp1 = dxkp1.^2 + dykp1.^2;
rki = 1./(rk + eye(Nvorts)) - eye(Nvorts); % should be zero along diagonal.  test.

gfunvals = gfun( (rkp1.*rki).^2 );
gfunvals(isinf(gfunvals)) = 0.;

gfunri2 = rki.^2 .* gfunvals;

xcoord = (xkp1-xk)/dt + i2p * (ymid .* gfunri2)*gvals;
ycoord = (ykp1-yk)/dt - i2p * (xmid .* gfunri2)*gvals;

system = [xcoord; ycoord];

%disp(sum(xcoord.^2 + ycoord.^2))

end

function val = gfun(z)

    if z == 1.
        val = 1.;
    else
        val = log(z)./(z-1);
    end

end