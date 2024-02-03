function radial_density = density_check(xpos, zpos, r0, lscl, ndists)
    Nvorts = size(xpos,2);
    cfun = @(x,z) sqrt(x.^2 + z.^2);
    radial_density = zeros(ndists,1);
    lscl = lscl/sqrt(8*pi);
    for jj=0:ndists-1
        rout = (r0-jj*lscl);
        rin = (r0-(jj+1)*lscl);
        vout = (cfun(xpos, zpos)<=rout);
        vin = (cfun(xpos, zpos)>rin);
        ref_inds = vout.*vin;
        radial_density(jj+1) = sum(ref_inds)/(pi*(rout.^2-rin.^2));
    end
    radial_density = flip(radial_density, 1);
end