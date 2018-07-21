function error_compare(mu,gam,omega,tf)

[tvals1,evec1] = waves_over_vortices_gen_curve(200,mu,gam,omega,tf,6);
[tvals2,evec2] = waves_over_vortices_gen_curve(300,mu,gam,omega,tf,6);
[tvals3,evec3] = waves_over_vortices_gen_curve(400,mu,gam,omega,tf,6);

plot(tvals1,evec1,'k-',tvals2,evec2,'k--',tvals3,evec3,'k-.','LineWidth',2)
h = set(gca,'FontSize',30);
set(h,'Interpreter','LaTeX')
xlabel('$t$','Interpreter','LaTeX','FontSize',30)
ylabel('$\mathcal{E}(t)$','Interpreter','LaTeX','FontSize',30)
legend({'$h=.01$','$h=.067$','$h=.05$'},'Interpreter','LaTeX')
