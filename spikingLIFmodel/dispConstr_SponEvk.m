%Comparing Spon&Evoked data sets

load netParms %!!! assuming same Parms_data (and thus id1_, id2_, etc !!!

% -- alter id[1/2]_ob to capture less I's --, arrays dont capture all granule cells
%used to alter which OB cells select (more granule sim but array records less)
%pct_[]  is % of TOTAL, so if 50%, then use all M/T and same number of granule
pct_granOB=0.5; %pct_granOB must be [0.02,0.80]
pct_mtOB=1-pct_granOB;
%indices of OB granule cells in correl(OB) calc
obGran_ind=sort(randperm(round(.8*Nob),round(pct_granOB/pct_mtOB*.2*Nob)))';
%indices of all OB (gran or M/T) in correl(OB) calc; use all M/T
ob_svInd=[(1:.2*Nob)'; obGran_ind+(.2*Nob)]; %first .2*Nob are excit, last .8*Nob are inhib granule
id1OBtrim=[]; 
id2OBtrim=[];
cnt=1;
for j=1:szRob
    if( sum((id1_ob(j)+1-ob_svInd)==0)>0 &&  sum((id2_ob(j)+1-ob_svInd)==0)>0) %need BOTH id1/id2 to have ob_svInd
        id1OBtrim=[id1OBtrim; id1_ob(j)+1]; %tack on
        id2OBtrim=[id2OBtrim; id2_ob(j)+1];
        cnt=cnt+1; %update count
    end
end


np_indx=(1:Npc)';
no_indx=(1:Nob)';

% Loading sims in Spontaneous state (black)
load d_spnHet
%change var names to compare
nuPC1=nuPC;
nuOB1=nuOB;
mn_PC1=mn_PC;
var_PC1=var_PC;
mn_OB1=mn_OB;
var_OB1=var_OB;
icov_pp1=icov_pp;
icov_oo1=icov_oo;

% Loading sims in Evoked state (red)
load d_evkHet
%change var names to compare
nuPC2=nuPC;
nuOB2=nuOB;
mn_PC2=mn_PC;
var_PC2=var_PC;
mn_OB2=mn_OB;
var_OB2=var_OB;
icov_pp2=icov_pp;
icov_oo2=icov_oo;

%Fano Factor of spike counts as a function of window size
FF_p1=var_PC1./mn_PC1;
FF_p1(isnan(FF_p1))=0; %pad 0 mn_E with 0's
FF_p2=var_PC2./mn_PC2;
FF_p2(isnan(FF_p2))=0; %pad 0 mn_E with 0's
FF_o1=var_OB1./mn_OB1;
FF_o1(isnan(FF_o1))=0; %pad 0 mn_I with 0's
FF_o2=var_OB2./mn_OB2;
FF_o2(isnan(FF_o2))=0; %pad 0 mn_I with 0's

%-- correlation/cov ----
rhoPP1=icov_pp1./(sqrt(var_PC1(id1_pc+1,:).*var_PC1(id2_pc+1,:)));
rhoPP1(isnan(rhoPP1))=0; %pad 0 var with 0's
rhoPP2=icov_pp2./(sqrt(var_PC2(id1_pc+1,:).*var_PC2(id2_pc+1,:)));
rhoPP2(isnan(rhoPP2))=0; %pad 0 var with 0's
%rhoOO1=icov_oo1./sqrt(var_OB1(id1OBtrim,:).*var_OB1(id2OBtrim,:));
rhoOO1=icov_oo1./sqrt(var_OB1(id1_ob+1,:).*var_OB1(id2_ob+1,:));
rhoOO1(isnan(rhoOO1))=0; %pad 0 var with 0's
%rhoOO2=icov_oo2./sqrt(var_OB2(id1OBtrim,:).*var_OB2(id2OBtrim,:));
rhoOO2=icov_oo2./sqrt(var_OB2(id1_ob+1,:).*var_OB2(id2_ob+1,:));
rhoOO2(isnan(rhoOO2))=0; %pad 0 var with 0's

%----Spontaneous Constraints ----

spon_PC_rt_less_OB=mean(nuPC1)<mean(nuOB1)

figure
hold on
plot(T_win./1000,mean(rhoPP1),'color',[0 .5 0],'LineWidth',4)
plot(T_win./1000,mean(rhoOO1),'b','LineWidth',4)
set(gca,'FontSize',18)
xlabel('Time Window (s)')
ylabel('Correlation')
title('Spontaneous State')
legend('PC','OB')

figure
hold on
plot(T_win./1000,mean(FF_p1),'color',[0 .5 0],'LineWidth',4)
plot(T_win./1000,mean(FF_o1),'b','LineWidth',4)
set(gca,'FontSize',18)
xlabel('Time Window (s)')
ylabel('Fano Factor')
title('Spontaneous State')
legend('PC','OB')

%----Evoked Constraints ----

evk_PC_rt_less_OB=mean(nuPC2)<mean(nuOB2)

figure
hold on
plot(T_win./1000,mean(rhoPP2),'color',[0 .5 0],'LineWidth',4)
plot(T_win./1000,mean(rhoOO2),'b','LineWidth',4)
set(gca,'FontSize',18)
xlabel('Time Window (s)')
ylabel('Correlation')
title('Evoked State')
legend('PC','OB')

figure
hold on
plot(T_win./1000,mean(icov_pp2)./(T_win'./1000),'color',[0 .5 0],'LineWidth',4)
plot(T_win./1000,mean(icov_oo2)./(T_win'./1000),'b','LineWidth',4)
set(gca,'FontSize',18)
xlabel('Time Window (s)')
ylabel('Covariance/T')
title('Evoked State')
legend('PC','OB')

figure
hold on
plot(T_win./1000,mean(var_PC2)./(T_win'./1000),'color',[0 0.5 0],'LineWidth',4)
plot(T_win./1000,mean(var_OB2)./(T_win'./1000),'b','LineWidth',4)
set(gca,'FontSize',18)
xlabel('Time Window (s)')
ylabel('Variance/T')
title('Evoked State')
legend('PC','OB')

%----Spont+Evk Constraints ----

OB_spLessEv=mean(nuOB1)<mean(nuOB2)


PC_spLessEv=mean(nuPC1)<mean(nuPC2)

figure
hold on
plot(T_win./1000,mean(var_OB1)./(T_win'./1000),'k','LineWidth',4)
plot(T_win./1000,mean(var_OB2)./(T_win'./1000),'r','LineWidth',4)
set(gca,'FontSize',18)
xlabel('Time Window (s)')
ylabel('Variance/T')
title('Olfactory Bulb')
legend('Spontaneous','Evoked')

figure
hold on
plot(T_win./1000,mean(FF_p1),'k','LineWidth',4)
plot(T_win./1000,mean(FF_p2),'r','LineWidth',4)
set(gca,'FontSize',18)
xlabel('Time Window (s)')
ylabel('Fano Factor')
title('Piriform Cortex')
legend('Spontaneous','Evoked')


figure
hold on
plot(T_win./1000,mean(rhoPP1),'k','LineWidth',4)
plot(T_win./1000,mean(rhoPP2),'r','LineWidth',4)
set(gca,'FontSize',18)
xlabel('Time Window (s)')
ylabel('Correlation')
title('Piriform Cortex')
legend('Spontaneous','Evoked')

mean(nuOB1)
std(nuOB1)
mean(nuPC1)
std(nuPC1)
mean(nuOB2)
std(nuOB2)
mean(nuPC2)
std(nuPC2)