%script to display stats as a function of Twin, COMBINING both experiments
%Very similar to dispCorTwin
%similar to dispTwinBoth but now scaling std by sqrt(N) so it is standard
%error

N_ob=41; %from both data sets
N_pc=73;
Npairs_ob=23*(23-1)*.5+18*(18-1)*.5;
Npairs_pc=38*(38-1)*.5+35*(35-1)*.5;
Npairs_PB=23*38+18*35;

load CoVar_Both

% PC: cf Spont and Evoked
figure(61) 
hold on
shadedErrorBar(Twin_vec,mnCrPC_sponBoth,stdCrPC_sponBoth./sqrt(Npairs_pc),'k',1)
shadedErrorBar(Twin_vec,mnCrPC_odBoth,stdCrPC_odBoth./sqrt(Npairs_pc),'r',1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'Spontanous=black, Evoked=red';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Piriform Cortex')
xlabel('Time Window (s)')
ylabel('Correlation')

% OB: cf Spont and Evoked
figure 
hold on
shadedErrorBar(Twin_vec,mnCrOB_sponBoth,stdCrOB_sponBoth./sqrt(Npairs_ob),'k',1)
shadedErrorBar(Twin_vec,mnCrOB_odBoth,stdCrOB_odBoth./sqrt(Npairs_ob),'r',1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'Spontanous=black, Evoked=red';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Olfactory Bulb')
xlabel('Time Window (s)')
ylabel('Correlation')

% Spon: cf OB and PC
figure(63) 
hold on
shadedErrorBar(Twin_vec,mnCrOB_sponBoth,stdCrOB_sponBoth./sqrt(Npairs_ob),'b',1)
shadedErrorBar(Twin_vec,mnCrPC_sponBoth,stdCrPC_sponBoth./sqrt(Npairs_pc),{'color',[0 .5 0]},1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'OB=blue, PC=green';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Spontaneous State')
xlabel('Time Window (s)')
ylabel('Correlation')

% Evoked: cf OB and PC
figure(64) 
hold on
shadedErrorBar(Twin_vec,mnCrOB_odBoth,stdCrOB_odBoth./sqrt(Npairs_ob),'b',1)
shadedErrorBar(Twin_vec,mnCrPC_odBoth,stdCrPC_odBoth./sqrt(Npairs_pc),{'color',[0 .5 0]},1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'OB=blue, PC=green';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Evoked State')
xlabel('Time Window (s)')
ylabel('Correlation')

% PC-OB: cf Spont and Evoked
figure 
hold on
shadedErrorBar(Twin_vec,mnCrPB_sponBoth,stdCrPB_sponBoth./sqrt(Npairs_PB),'m',1)
shadedErrorBar(Twin_vec,mnCrPB_odBoth,stdCrPB_odBoth./sqrt(Npairs_PB),'c',1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'Spontanous=magenta, Evoked=cyan';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
xlabel('Time Window (s)')
ylabel('PC-OB Corrleation')


%Fano Factor of individual cells
% PC: cf Spont and Evoked
figure(62)
hold on
shadedErrorBar(Twin_vec,mnFFpc_sponBoth,stdFFpc_sponBoth./sqrt(N_pc),'k',1)
shadedErrorBar(Twin_vec,mnFFpc_odBoth,stdFFpc_odBoth./sqrt(N_pc),'r',1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'Spontanous=black, Evoked=red';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Piriform Cortex')
xlabel('Time Window (s)')
ylabel('Fano Factor')

% OB: cf Spont and Evoked
figure
hold on
shadedErrorBar(Twin_vec,mnFFob_sponBoth,stdFFob_sponBoth./sqrt(N_ob),'k',1)
shadedErrorBar(Twin_vec,mnFFob_odBoth,stdFFob_odBoth./sqrt(N_ob),'r',1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'Spontanous=black, Evoked=red';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Olfactory Bulb')
xlabel('Time Window (s)')
ylabel('Fano Factor')

% Spont: cf PC and OB
figure(65)
hold on
shadedErrorBar(Twin_vec,mnFFob_sponBoth,stdFFob_sponBoth./sqrt(N_ob),'b',1)
shadedErrorBar(Twin_vec,mnFFpc_sponBoth,stdFFpc_sponBoth./sqrt(N_pc),{'color',[0 .5 0]},1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'OB=blue, PC=green';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Spontaneous State')
xlabel('Time Window (s)')
ylabel('Fano Factor')

% Evoked: cf PC and OB
figure
hold on
shadedErrorBar(Twin_vec,mnFFob_odBoth,stdFFob_odBoth./sqrt(N_ob),'b',1)
shadedErrorBar(Twin_vec,mnFFpc_odBoth,stdFFpc_odBoth./sqrt(N_pc),{'color',[0 .5 0]},1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'OB=blue, PC=green';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Evoked State')
xlabel('Time Window (s)')
ylabel('Fano Factor')

%________Just checking the Vars________
figure
hold on
shadedErrorBar(Twin_vec,mnVrpc_sponBoth./Twin_vec,stdVrpc_sponBoth./sqrt(N_pc)./Twin_vec,'k',1)
shadedErrorBar(Twin_vec,mnVrpc_odBoth./Twin_vec,stdVrpc_odBoth./sqrt(N_pc)./Twin_vec,'r',1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'Spontanous=black, Evoked=red';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Piriform Cortex')
xlabel('Time Window (s)')
ylabel('Variance/T')

% OB: cf Spont and Evoked
figure(68)
hold on
shadedErrorBar(Twin_vec,mnVrob_sponBoth./Twin_vec,stdVrob_sponBoth./sqrt(N_ob)./Twin_vec,'k',1)
shadedErrorBar(Twin_vec,mnVrob_odBoth./Twin_vec,stdVrob_odBoth./sqrt(N_ob)./Twin_vec,'r',1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'Spontanous=black, Evoked=red';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Olfactory Bulb')
xlabel('Time Window (s)')
ylabel('Variance/T')

% Spont: cf PC and OB
figure
hold on
shadedErrorBar(Twin_vec,mnVrob_sponBoth./Twin_vec,stdVrob_sponBoth./sqrt(N_ob)./Twin_vec,'b',1)
shadedErrorBar(Twin_vec,mnVrpc_sponBoth./Twin_vec,stdVrpc_sponBoth./sqrt(N_pc)./Twin_vec,{'color',[0 .5 0]},1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'OB=blue, PC=green';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Spontaneous State')
xlabel('Time Window (s)')
ylabel('Variance/T')

% Evoked: cf PC and OB
figure(66)
hold on
shadedErrorBar(Twin_vec,mnVrob_odBoth./Twin_vec,stdVrob_odBoth./sqrt(N_ob)./Twin_vec,'b',1)
shadedErrorBar(Twin_vec,mnVrpc_odBoth./Twin_vec,stdVrpc_odBoth./sqrt(N_pc)./Twin_vec,{'color',[0 .5 0]},1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'OB=blue, PC=green';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Evoked State')
xlabel('Time Window (s)')
ylabel('Variance/T')


%check Cov for same trends

% PC: cf Spont and Evoked
figure
hold on
shadedErrorBar(Twin_vec,mnCovPC_sponBoth./Twin_vec,stdCovPC_sponBoth./sqrt(Npairs_pc)./Twin_vec,'k',1)
shadedErrorBar(Twin_vec,mnCovPC_odBoth./Twin_vec,stdCovPC_odBoth./sqrt(Npairs_pc)./Twin_vec,'r',1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'Spontanous=black, Evoked=red';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Piriform Cortex')
xlabel('Time Window (s)')
ylabel('Covariance/T')

% OB: cf Spont and Evoked
figure
hold on
shadedErrorBar(Twin_vec,mnCovOB_sponBoth./Twin_vec,stdCovOB_sponBoth./sqrt(Npairs_ob)./Twin_vec,'k',1)
shadedErrorBar(Twin_vec,mnCovOB_odBoth./Twin_vec,stdCovOB_odBoth./sqrt(Npairs_ob)./Twin_vec,'r',1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'Spontanous=black, Evoked=red';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Olfactory Bulb')
xlabel('Time Window (s)')
ylabel('Covariance/T')

% Spon: cf OB and PC
figure
hold on
shadedErrorBar(Twin_vec,mnCovOB_sponBoth./Twin_vec,stdCovOB_sponBoth./sqrt(Npairs_ob)./Twin_vec,'b',1)
shadedErrorBar(Twin_vec,mnCovPC_sponBoth./Twin_vec,stdCovPC_sponBoth./sqrt(Npairs_pc)./Twin_vec,{'color',[0 .5 0]},1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'OB=blue, PC=green';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Spontaneous State')
xlabel('Time Window (s)')
ylabel('Covariance/T')

% Evoked: cf OB and PC
figure(67)
hold on
shadedErrorBar(Twin_vec,mnCovOB_odBoth./Twin_vec,stdCovOB_odBoth./sqrt(Npairs_ob)./Twin_vec,'b',1)
shadedErrorBar(Twin_vec,mnCovPC_odBoth./Twin_vec,stdCovPC_odBoth./sqrt(Npairs_pc)./Twin_vec,{'color',[0 .5 0]},1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'OB=blue, PC=green';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
title('Evoked State')
xlabel('Time Window (s)')
ylabel('Covariance/T')

% PC-OB: cf Spont and Evoked
figure
hold on
shadedErrorBar(Twin_vec,mnCovPB_sponBoth./Twin_vec,stdCovPB_sponBoth./sqrt(Npairs_PB)./Twin_vec,'m',1)
shadedErrorBar(Twin_vec,mnCovPB_odBoth./Twin_vec,stdCovPB_odBoth./sqrt(Npairs_PB)./Twin_vec,'c',1)
set(gca,'FontSize',18)
dim = [.5 .6 .3 .3];
str1 = 'Spontanous=magenta, Evoked=cyan';
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
xlabel('Time Window (s)')
ylabel('PC-OB Covariance')

