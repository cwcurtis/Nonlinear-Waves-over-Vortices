function res=E2(z,m,gam,ep)

    terf = erf(abs(ep*pi*m/2-gam*z/ep));
    inds = z < ep^2*pi*m/(2*gam);
    inds = 2*inds - 1;
    res = 1 + inds.*terf;

end