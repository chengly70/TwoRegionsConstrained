%commands to plot the stats BY ODOR; 
%plotting the Corr & Cov of all 3 pairs (PC,OB,PC-OB) in Spont. vs Evoked

which_set=input('Enter 0 (exp1+2) or 1 (exp3+4): '); 
if(which_set==0)
    load CoVar_Exp1+2
else
    load CoVar_Exp3+4
end

figure
hold on
plot(Twin_vec,mnCrPC_sponBoth,'k','LineWidth',4)
plot(Twin_vec,mnCrPC_spon1,'r','LineWidth',4)
plot(Twin_vec,mnCrPC_spon2,'b','LineWidth',4)
plot(Twin_vec,mnCrPC_odBoth,'k','LineWidth',4)
plot(Twin_vec,mnCrPC_od1,'r','LineWidth',4)
plot(Twin_vec,mnCrPC_od2,'b','LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('PC Correlation')

figure
hold on
plot(Twin_vec,mnCrOB_sponBoth,'k','LineWidth',4)
plot(Twin_vec,mnCrOB_spon1,'r','LineWidth',4)
plot(Twin_vec,mnCrOB_spon2,'b','LineWidth',4)
plot(Twin_vec,mnCrOB_odBoth,'k','LineWidth',4)
plot(Twin_vec,mnCrOB_od1,'r','LineWidth',4)
plot(Twin_vec,mnCrOB_od2,'b','LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('OB Correlation')

figure
hold on
plot(Twin_vec,mnCrPB_sponBoth,'k','LineWidth',4)
plot(Twin_vec,mnCrPB_spon1,'r','LineWidth',4)
plot(Twin_vec,mnCrPB_spon2,'b','LineWidth',4)
plot(Twin_vec,mnCrPB_odBoth,'k','LineWidth',4)
plot(Twin_vec,mnCrPB_od1,'r','LineWidth',4)
plot(Twin_vec,mnCrPB_od2,'b','LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('PC-OB Correlation')

figure
hold on
plot(Twin_vec,mnCovPC_sponBoth./Twin_vec,'k','LineWidth',4)
plot(Twin_vec,mnCovPC_spon1./Twin_vec,'r','LineWidth',4)
plot(Twin_vec,mnCovPC_spon2./Twin_vec,'b','LineWidth',4)
plot(Twin_vec,mnCovPC_odBoth./Twin_vec,'k','LineWidth',4)
plot(Twin_vec,mnCovPC_od1./Twin_vec,'r','LineWidth',4)
plot(Twin_vec,mnCovPC_od2./Twin_vec,'b','LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('PC Covariance/T')

figure
hold on
plot(Twin_vec,mnCovOB_sponBoth./Twin_vec,'k','LineWidth',4)
plot(Twin_vec,mnCovOB_spon1./Twin_vec,'r','LineWidth',4)
plot(Twin_vec,mnCovOB_spon2./Twin_vec,'b','LineWidth',4)
plot(Twin_vec,mnCovOB_odBoth./Twin_vec,'k','LineWidth',4)
plot(Twin_vec,mnCovOB_od1./Twin_vec,'r','LineWidth',4)
plot(Twin_vec,mnCovOB_od2./Twin_vec,'b','LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('OB Covariance/T')

figure
hold on
plot(Twin_vec,mnCovPB_sponBoth./Twin_vec,'k','LineWidth',4)
plot(Twin_vec,mnCovPB_spon1./Twin_vec,'r','LineWidth',4)
plot(Twin_vec,mnCovPB_spon2./Twin_vec,'b','LineWidth',4)
plot(Twin_vec,mnCovPB_odBoth./Twin_vec,'k','LineWidth',4)
plot(Twin_vec,mnCovPB_od1./Twin_vec,'r','LineWidth',4)
plot(Twin_vec,mnCovPB_od2./Twin_vec,'b','LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('PC-OB Covariance/T')