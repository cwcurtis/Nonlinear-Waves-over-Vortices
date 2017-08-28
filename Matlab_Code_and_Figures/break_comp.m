function break_comp(tf)

K = 256;

fvals = linspace(.01,.3,32);
tvalspt1 = zeros(length(fvals),1);
tvalspt2 = zeros(length(fvals),1);
fvalspt1 = zeros(length(fvals),1);
fvalspt2 = zeros(length(fvals),1);

parfor jj=1:length(fvals)
    [tvalspt1(jj),fvalspt1(jj)] = waves_over_vortices_break(K,.1,sqrt(.1),fvals(jj),tf);
    [tvalspt2(jj),fvalspt2(jj)] = waves_over_vortices_break(K,.2,sqrt(.2),fvals(jj),tf);
end

clf

figure(1)
plot(fvals,tvalspt1,'k--',fvals,tvalspt2,'k','LineWidth',2)
set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
xlabel('F','FontName','Helvetica','FontSize',30,'FontWeight','bold')
ylabel('t_{b}','FontName','Helvetica','FontSize',30,'FontWeight','bold')

figure(2)
plot(fvals,fvalspt1,'k--',fvals,fvalspt2,'k','LineWidth',2)
set(gca,'FontSize',30,'FontName','Helvetica','FontWeight','bold')
xlabel('F','FontName','Helvetica','FontSize',30,'FontWeight','bold')
ylabel('\delta_{b}','FontName','Helvetica','FontSize',30,'FontWeight','bold')
