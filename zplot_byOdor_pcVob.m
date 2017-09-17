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
plot(Twin_vec,mnCrPC_sponBoth,'color',[0 0.5 0],'LineWidth',4)
plot(Twin_vec,mnCrPC_spon1,'color',[0 0.35 0],'LineWidth',4)
plot(Twin_vec,mnCrPC_spon2,'color',[0 0.7 0],'LineWidth',4)
plot(Twin_vec,mnCrOB_sponBoth,'color',[0 0 1],'LineWidth',4)
plot(Twin_vec,mnCrOB_spon1,'color',[0 0 0.7],'LineWidth',4)
plot(Twin_vec,mnCrOB_spon2,'color',[0 0 0.35],'LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('Spontaneous Correlation')

figure
hold on
plot(Twin_vec,mnCrPC_odBoth,'color',[0 0.5 0],'LineWidth',4)
plot(Twin_vec,mnCrPC_od1,'color',[0 0.35 0],'LineWidth',4)
plot(Twin_vec,mnCrPC_od2,'color',[0 0.7 0],'LineWidth',4)
plot(Twin_vec,mnCrOB_odBoth,'color',[0 0 1],'LineWidth',4)
plot(Twin_vec,mnCrOB_od1,'color',[0 0 0.7],'LineWidth',4)
plot(Twin_vec,mnCrOB_od2,'color',[0 0 0.35],'LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('Evoked Correlation')


figure
hold on
plot(Twin_vec,mnCovPC_sponBoth./Twin_vec,'color',[0 0.5 0],'LineWidth',4)
plot(Twin_vec,mnCovPC_spon1./Twin_vec,'color',[0 0.35 0],'LineWidth',4)
plot(Twin_vec,mnCovPC_spon2./Twin_vec,'color',[0 0.7 0],'LineWidth',4)
plot(Twin_vec,mnCovOB_sponBoth./Twin_vec,'color',[0 0 1],'LineWidth',4)
plot(Twin_vec,mnCovOB_spon1./Twin_vec,'color',[0 0 0.7],'LineWidth',4)
plot(Twin_vec,mnCovOB_spon2./Twin_vec,'color',[0 0 0.35],'LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('Spontaneous Covariance/T')

figure
hold on
plot(Twin_vec,mnCovPC_odBoth./Twin_vec,'color',[0 0.5 0],'LineWidth',4)
plot(Twin_vec,mnCovPC_od1./Twin_vec,'color',[0 0.35 0],'LineWidth',4)
plot(Twin_vec,mnCovPC_od2./Twin_vec,'color',[0 0.7 0],'LineWidth',4)
plot(Twin_vec,mnCovOB_odBoth./Twin_vec,'color',[0 0 1],'LineWidth',4)
plot(Twin_vec,mnCovOB_od1./Twin_vec,'color',[0 0 0.7],'LineWidth',4)
plot(Twin_vec,mnCovOB_od2./Twin_vec,'color',[0 0 0.35],'LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('Evoked Covariance/T')