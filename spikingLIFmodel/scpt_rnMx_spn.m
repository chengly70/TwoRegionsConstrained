%m-file to run the mex function, mx_jCorv_LIF2Pop_delay.c for Spont state
%a Monte Carlo sim to get all stats (frate, var, rndSmpl corr/cov) and error bars

Nrlz=50000;

%command to load inputs, connect, Parameters, etc
%load inputs, connect, Parameters, etc
load netParms
cr_pi=[0.5; 0.8]; %[crOB;crPC]
sigp=[0.05;0.1]; %[sigOB;sigPC]
Iap_OB=0.6*ones(size(ThresOB)); %Spontaneous
Iap_PC=0*ones(size(ThresPC));   %Spontaneous
%Iap_OB=0.9*ones(size(ThresOB));  %Evoked
%Iap_PC=0.4*ones(size(ThresPC));  %Evoked

% --- specify conductance strengths ----
gIO=7;   %I->E in OB
gieo=4;    %E->I in OB
geeo=2;   %E->E in OB
giio=2; %6;    %I->I in OB
gIP=20;   %I->E in PC
giep=8;  %4;  %E->I in PC
geep=5;  %2;  %E->E in PC
giip=6;    %I->I in PC
gCrsoep=0;%1; %cross; PC to OB(E)
gEP=15;     %cross; PC to OB(I)
gCrspio=0; %1; %cross; OB to PC(I)
gEO=10;     %cross; OB to PC(E)
gRaw_vec=[gIO;gieo;geeo;giio;gIP;giep;geep;giip;gCrsoep;gEP;gEO;gCrspio];
%normalize syn-strength by fraction of connected
g_IO=gIO/(Nob*fr_ItoEob);
g_ieo=gieo/(Nob*fr_EtoIob);
g_eeo=geeo/(Nob*fr_EtoEob);
g_iio=giio/(Nob*fr_ItoIob);
g_IP=gIP/(Npc*fr_ItoEpc);
g_iep=giep/(Npc*fr_EtoIpc);
g_eep=geep/(Npc*fr_EtoEpc);
g_iip=giip/(Npc*fr_ItoIpc);
g_Crsoep=gCrsoep/(Npc*fr_pcToE);
g_EP=gEP/(Npc*fr_pcToI);
g_Crspio=gCrspio/(Nob*fr_obToI);
g_EO=gEO/(Nob*fr_obToE);


g_vec=[g_IO;g_ieo;g_eeo;g_iio;g_IP;g_iep;g_eep;g_iip;g_Crsoep;g_EP;g_EO;g_Crspio]; %!!must be a 12x1 vector!!

% check that numbers won't crash C/mex 
if(mod(length(ThresOB)*0.2,1) ||  mod(length(ThresPC)*0.2,1) || length(id1_ob)~=length(id2_ob) || length(id1_pc)~=length(id2_pc))
    disp('Check parms, errors!');
    return;
end
if(size(W_oo,1)~=size(W_op,1) || size(W_pp,1)~=size(W_po,1) )
    disp('Coupling matrices are wrong size!');
    return;
end

%must match c-file mx_jCorv_LIF2Pop_delay.c!
T_win=[50; 100; 200; 250; 400; 500; 1000; 2000];

tic
    [nuPC,nuOB,mn_PC,mn_OB,var_PC,var_OB,icov_pp,icov_oo]=mx_jCorv_LIF2Pop_delay(W_oo,W_pp,W_op,W_po,g_vec,id1_ob,id2_ob,id1_pc,id2_pc,ThresOB,ThresPC,cr_pi,sigp,Iap_OB,Iap_PC,Nrlz);
toc


save d_spnHet nuPC nuOB mn_* var_* icov_* T_win sigp cr_pi Iap_*