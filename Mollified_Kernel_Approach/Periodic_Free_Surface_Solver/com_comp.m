function xcm = com_comp(ep,gam,xpos,zpos,gvals)

    chi = @(r) 1./pi.*(2.*exp(-(r.^2))-.5*exp(-(r/sqrt(2)).^2));
    diffmat = sqrt((xpos*ones(1,length(xpos))-ones(length(xpos),1)*xpos').^2 + gam^2*(zpos*ones(1,length(zpos))-ones(length(zpos),1)*zpos').^2);
    ivec = gam*chi(diffmat/ep)*gvals/(ep^2);
    Finterp = scatteredInterpolant(xpos,zpos,ivec);
    xmin = min(xpos);
    xmax = max(xpos);
    ddx = (xmax-xmin)/100;
    zmin = min(zpos);
    zmax = max(zpos);
    ddz = (zmax-zmin)/100;
    [Xmm,Zmm] = meshgrid((xmin:ddx:xmax),(zmin:ddz:zmax));
    Fvals = Finterp(Xmm,Zmm);
    intgrnd1 = Xmm.*Fvals;
    intgrnd2 = Fvals;
    xcm = ( intgrnd1(1,1) + intgrnd1(1,end) + intgrnd1(end,1) + intgrnd1(end,end) + 2*sum(intgrnd1(1,2:end-1)' + intgrnd1(end,2:end-1)' + intgrnd1(2:end-1,1) + intgrnd1(2:end-1,end)) + 4*sum(sum(intgrnd1(2:end-1,2:end-1))) )/...
          ( intgrnd2(1,1) + intgrnd2(1,end) + intgrnd2(end,1) + intgrnd2(end,end) + 2*sum(intgrnd2(1,2:end-1)' + intgrnd2(end,2:end-1)' + intgrnd2(2:end-1,1) + intgrnd2(2:end-1,end)) + 4*sum(sum(intgrnd2(2:end-1,2:end-1))) );