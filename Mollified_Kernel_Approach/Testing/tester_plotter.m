function tester_plotter(Nvorts,gam,ep,Ntrunc)

Nvals = 2.^(1:Ntrunc);

Kfnorm = zeros(Ntrunc-1,1);
Ksnorm = zeros(Ntrunc-1,1);

Kmatx_fp = kernel_tester(Nvorts,gam,ep,Nvals(1));
Kmatx_sp = kernel_tester_sum(Nvorts,gam,ep,Nvals(1));

for jj = 2:Ntrunc
   Kmatx_f = kernel_tester(Nvorts,gam,ep,Nvals(jj));
   Kfnorm(jj-1) = norm(Kmatx_f - Kmatx_fp);
   Kmatx_fp = Kmatx_f;
end

for jj = 2:Ntrunc
   Kmatx_s = kernel_tester_sum(Nvorts,gam,ep,Nvals(jj));
   Ksnorm(jj-1) = norm(Kmatx_s - Kmatx_sp);
   Kmatx_sp = Kmatx_s;
end

disp(num2str(abs(Kmatx_sp-Kmatx_fp)))

figure(1)
plot((2:Ntrunc)',log10(Kfnorm))

figure(2)
plot((2:Ntrunc)',log10(Ksnorm))

end

