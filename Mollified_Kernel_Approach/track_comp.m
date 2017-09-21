function track_comp(K,mu,gam,tf)

dt = 1e-2;

[xtracknv,ztracknv] = waves_over_vortices_tracks(K,mu,gam,0,tf,dt);
[xtrackwk,ztrackwk] = waves_over_vortices_tracks(K,mu,gam,.02,tf,dt);
[xtrackst,ztrackst] = waves_over_vortices_tracks(K,mu,gam,.2,tf,dt);

figure(1)
hold on 
plot(xtracknv(1,:),ztracknv(1,:),'k--',xtrackwk(1,:),ztrackwk(1,:),'k-.',xtrackst(1,:),ztrackst(1,:),'k','LineWidth',2)
plot(xtracknv(1,1),ztracknv(1,1),'.','MarkerSize',26','color',[0.8 0.8 0.8]);
plot(xtracknv(1,end),ztracknv(1,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
plot(xtrackwk(1,end),ztrackwk(1,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
plot(xtrackst(1,end),ztrackst(1,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
hold off  
set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
ylabel('z','FontName','Helvetica','FontSize',30,'FontWeight','bold')

figure(2)
hold on 
plot(xtracknv(2,:),ztracknv(2,:),'k--',xtrackwk(2,:),ztrackwk(2,:),'k-.',xtrackst(2,:),ztrackst(2,:),'k','LineWidth',2)
plot(xtracknv(2,1),ztracknv(2,1),'.','MarkerSize',26','color',[0.8 0.8 0.8]);
plot(xtracknv(2,end),ztracknv(2,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
plot(xtrackwk(2,end),ztrackwk(2,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
plot(xtrackst(2,end),ztrackst(2,end),'.','MarkerSize',26','color',[0.1 0.1 0.1]);
hold off  
set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
xlabel('x','FontName','Helvetica','FontSize',30,'FontWeight','bold')
ylabel('z','FontName','Helvetica','FontSize',30,'FontWeight','bold')
