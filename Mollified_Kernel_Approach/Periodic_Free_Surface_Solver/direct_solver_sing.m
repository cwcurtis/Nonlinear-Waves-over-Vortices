function nl = direct_solver_sing(zpos,gval,Mx)
p2M = pi/(2*Mx);
%dzp = 2*zpos;
%dzzp = 1i*dzp;
pkerp = cot(p2M*(2*1i*zpos));
Kval = -pkerp*gval;
nl = [imag(Kval) real(Kval)];