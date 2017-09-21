function conv_comp(mu,gam,tf)

Kvals = [32 64 128 256];
K = 256;
KT = 2*K;

[sval1,~,xpos1,zpos1] = waves_over_vortices_verify(K,mu,gam,.1,tf);
[sval2,~,xpos2,zpos2] = waves_over_vortices_verify(K,mu,gam,.2,tf);
[sval3,~,xpos3,zpos3] = waves_over_vortices_verify(K,mu,gam,.3,tf);

Xmeshfin = linspace(-1,1,KT+1);
Xmeshfin = Xmeshfin(1:KT)';

clf
hold on 

%plot(Xmesh,log10(abs(sval1-sval2(1:2:end))./abs(sval1)),'k','LineWidth',2)
plot(Xmeshfin,sval1,'k-.','LineWidth',2)
plot(Xmeshfin,sval2,'k--','LineWidth',2)
plot(Xmeshfin,sval3,'k','LineWidth',2)

hold off

set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
%ylabel('log_{10}|\eta_{512}-\eta_{256}|','FontName','Helvetica','FontSize',30,'FontWeight','bold')
