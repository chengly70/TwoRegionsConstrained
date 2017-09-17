%commands to plot the stats BY ODOR;
%plotting the FF & Var of 2 types in Spont. vs Evoked

which_set=input('Enter 0 (exp1+2) or 1 (exp3+4): '); 
if(which_set==0)
    load CoVar_Exp1+2
else
    load CoVar_Exp3+4
end

figure
hold on
plot(Twin_vec,mnFFpc_sponBoth,'k','LineWidth',4)
plot(Twin_vec,mnFFpc_spon1,'r','LineWidth',4)
plot(Twin_vec,mnFFpc_spon2,'b','LineWidth',4)
plot(Twin_vec,mnFFpc_odBoth,'k','LineWidth',4)
plot(Twin_vec,mnFFpc_od1,'r','LineWidth',4)
plot(Twin_vec,mnFFpc_od2,'b','LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('PC Fano Factor')

figure
hold on
plot(Twin_vec,mnFFob_sponBoth,'k','LineWidth',4)
plot(Twin_vec,mnFFob_spon1,'r','LineWidth',4)
plot(Twin_vec,mnFFob_spon2,'b','LineWidth',4)
plot(Twin_vec,mnFFob_odBoth,'k','LineWidth',4)
plot(Twin_vec,mnFFob_od1,'r','LineWidth',4)
plot(Twin_vec,mnFFob_od2,'b','LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('OB Fano Factor')


figure
hold on
plot(Twin_vec,mnVrpc_sponBoth./Twin_vec,'k','LineWidth',4)
plot(Twin_vec,mnVrpc_spon1./Twin_vec,'r','LineWidth',4)
plot(Twin_vec,mnVrpc_spon2./Twin_vec,'b','LineWidth',4)
plot(Twin_vec,mnVrpc_odBoth./Twin_vec,'k','LineWidth',4)
plot(Twin_vec,mnVrpc_od1./Twin_vec,'r','LineWidth',4)
plot(Twin_vec,mnVrpc_od2./Twin_vec,'b','LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('PC Variance/T')

figure
hold on
plot(Twin_vec,mnVrob_sponBoth./Twin_vec,'k','LineWidth',4)
plot(Twin_vec,mnVrob_spon1./Twin_vec,'r','LineWidth',4)
plot(Twin_vec,mnVrob_spon2./Twin_vec,'b','LineWidth',4)
plot(Twin_vec,mnVrob_odBoth./Twin_vec,'k','LineWidth',4)
plot(Twin_vec,mnVrob_od1./Twin_vec,'r','LineWidth',4)
plot(Twin_vec,mnVrob_od2./Twin_vec,'b','LineWidth',4)
box off
set(gca,'FontSize',18)
xlabel('Twin (s)')
ylabel('OB Variance/T')

