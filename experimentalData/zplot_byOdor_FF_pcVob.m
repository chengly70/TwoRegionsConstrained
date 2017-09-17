%commands to plot the stats BY ODOR;
%plotting the Corr & Cov of PC vs. OB pairs

which_set=input('Enter 0 (exp1+2) or 1 (exp3+4): '); 
if(which_set==0)
    load CoVar_Exp1+2
else
    load CoVar_Exp3+4
end

figure
hold on
plot(Twin_vec,mnFFpc_sponBoth,'color',[0 0.5 0],'LineWidth',4)
plot(Twin_vec,mnFFpc_spon1,'color',[0 0.35 0],'LineWidth',4)
plot(Twin_vec,mnFFpc_spon2,'color',[0 0.7 0],'LineWidth',4)
plot(Twin_vec,mnFFob_sponBoth,'color',[0 0 1],'LineWidth',4)
plot(Twin_vec,mnFFob_spon1,'color',[0 0 0.7],'LineWidth',4)
plot(Twin_vec,mnFFob_spon2,'color',[0 0 0.35],'LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('Spontaneous Fano Factor')

figure
hold on
plot(Twin_vec,mnFFpc_odBoth,'color',[0 0.5 0],'LineWidth',4)
plot(Twin_vec,mnFFpc_od1,'color',[0 0.35 0],'LineWidth',4)
plot(Twin_vec,mnFFpc_od2,'color',[0 0.7 0],'LineWidth',4)
plot(Twin_vec,mnFFob_odBoth,'color',[0 0 1],'LineWidth',4)
plot(Twin_vec,mnFFob_od1,'color',[0 0 0.7],'LineWidth',4)
plot(Twin_vec,mnFFob_od2,'color',[0 0 0.35],'LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('Evoked Fano Factor')


figure
hold on
plot(Twin_vec,mnVrpc_sponBoth./Twin_vec,'color',[0 0.5 0],'LineWidth',4)
plot(Twin_vec,mnVrpc_spon1./Twin_vec,'color',[0 0.35 0],'LineWidth',4)
plot(Twin_vec,mnVrpc_spon2./Twin_vec,'color',[0 0.7 0],'LineWidth',4)
plot(Twin_vec,mnVrob_sponBoth./Twin_vec,'color',[0 0 1],'LineWidth',4)
plot(Twin_vec,mnVrob_spon1./Twin_vec,'color',[0 0 0.7],'LineWidth',4)
plot(Twin_vec,mnVrob_spon2./Twin_vec,'color',[0 0 0.35],'LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('Spontaneous Variance/T')

figure
hold on
plot(Twin_vec,mnVrpc_odBoth./Twin_vec,'color',[0 0.5 0],'LineWidth',4)
plot(Twin_vec,mnVrpc_od1./Twin_vec,'color',[0 0.35 0],'LineWidth',4)
plot(Twin_vec,mnVrpc_od2./Twin_vec,'color',[0 0.7 0],'LineWidth',4)
plot(Twin_vec,mnVrob_odBoth./Twin_vec,'color',[0 0 1],'LineWidth',4)
plot(Twin_vec,mnVrob_od1./Twin_vec,'color',[0 0 0.7],'LineWidth',4)
plot(Twin_vec,mnVrob_od2./Twin_vec,'color',[0 0 0.35],'LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('Evoked Variance/T')