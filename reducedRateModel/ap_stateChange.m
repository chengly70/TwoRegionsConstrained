%Checking constraints across states (Sp --> Ev)

% ___ with crPC=0.35<crOB=0.3 ___
load dSp_st

% changing the names of the stats in Spont state (first) so can compare
frCorr_obSpon=obCorrRt_s;
frCov_obSpon=obCovRt_s;
frCorr_pcSpon=pcCorrRt_s;
frCov_pcSpon=pcCovRt_s;
frMn_Spon=frate_s; %all 6
frVr_Spon=rtVar_s; %all 6

%assuming same lengths and same d_g
d_g=gEP(2)-gEP(1); 
len_g=length(gEP);
leng_2=len_g^2;
leng_3=len_g^3;
tot_gs=length(gs_valid);

%calc mean stats (across 3, small sample)
mn_ob_frSpon=mean(frMn_Spon(:,1:3),2);
mn_pc_frSpon=mean(frMn_Spon(:,4:6),2);
vr_ob_frSpon=mean(frVr_Spon(:,1:3),2);
vr_pc_frSpon=mean(frVr_Spon(:,4:6),2);
mnFano_pc_Spon=mean(frVr_Spon(:,4:6)./frMn_Spon(:,4:6),2);
mnFano_ob_Spon=mean(frVr_Spon(:,1:3)./frMn_Spon(:,1:3),2);
%!!! in next 4 variables, don't count params that have gs_valid==1 (& others!!)
mnCorr_ob_frSpon=mean(frCorr_obSpon,2); 
mnCorr_pc_frSpon=mean(frCorr_pcSpon,2);
mnCov_ob_frSpon=mean(frCov_obSpon,2);
mnCov_pc_frSpon=mean(frCov_pcSpon,2);

feas_reg=(gs_valid>=3); %logical vector indicating feasible region (1) vs infeas (0)
%--- Constraint: rate PC < rate OB in Spont  ------
spConst_nu_pcVob = ((mn_pc_frSpon<mn_ob_frSpon).*feas_reg)==1;
%--- Constraint: FF of PC > OB in Spont  ------
spConst_fano_pcVob = ((mnFano_pc_Spon > mnFano_ob_Spon).*feas_reg)==1; 

%--- Constraint: Correl of PC > OB in Spont  ------
spConst_corr_pcVob = ((mnCorr_pc_frSpon > mnCorr_ob_frSpon).*feas_reg)==1;

feas_Spont=(spConst_fano_pcVob.*spConst_corr_pcVob.*spConst_nu_pcVob)==1;


pctSp_R_pcVob=sum(feas_reg)/tot_gs
pctSp_fano_pcVob=sum(spConst_fano_pcVob)/tot_gs
pctSp_corr_pcVob=sum(spConst_corr_pcVob)/tot_gs
pctSp_all=sum(feas_Spont)/tot_gs
pause


%------ Repeat above but now for evoked state --------

% ___ with crPC=0.35<crOB=0.3 ___
load dEv_st


% changing the names of the stats in Evoked state so can compare
frCorr_obEvk=obCorrRt_s;
frCov_obEvk=obCovRt_s;
frCorr_pcEvk=pcCorrRt_s;
frCov_pcEvk=pcCovRt_s;
frMn_Evk=frate_s; %all 6
frVr_Evk=rtVar_s; %all 6

%calc mean stats (across 3, small sample)
mn_ob_frEvk=mean(frMn_Evk(:,1:3),2);
mn_pc_frEvk=mean(frMn_Evk(:,4:6),2);
vr_ob_frEvk=mean(frVr_Evk(:,1:3),2);
vr_pc_frEvk=mean(frVr_Evk(:,4:6),2);
mnFano_pc_Evk=mean(frVr_Evk(:,4:6)./frMn_Evk(:,4:6),2);
mnFano_ob_Evk=mean(frVr_Evk(:,1:3)./frMn_Evk(:,1:3),2);
%!!! in next 4 variables, don't count params that have gs_valid==1 (& others!!)
mnCorr_ob_frEvk=mean(frCorr_obEvk,2); 
mnCorr_pc_frEvk=mean(frCorr_pcEvk,2);
mnCov_ob_frEvk=mean(frCov_obEvk,2);
mnCov_pc_frEvk=mean(frCov_pcEvk,2);

feas_reg=(gs_valid>=3); %logical vector indicating feasible region (1) vs infeas (0)
%--- Constraint: rate PC < rate OB in Evoked  ------
evConst_nu_pcVob = ((mn_pc_frEvk<mn_ob_frEvk).*feas_reg)==1;
%--- Constraint: Cov of PC < OB in Evoked  ------
evConst_Cov_pcVob = ((mnCov_pc_frEvk < mnCov_ob_frEvk).*feas_reg)==1;

%--- Constraint: Correl of PC < OB in Evoked  ------
evConst_corr_pcVob = ((mnCorr_pc_frEvk < mnCorr_ob_frEvk).*feas_reg)==1;

%--- Constraint: Var PC < OB in Evoked  ------
evConst_Vr_pcVob = ((vr_pc_frEvk < vr_ob_frEvk).*feas_reg)==1;

feas_Evk=(evConst_corr_pcVob.*evConst_Cov_pcVob.*evConst_Vr_pcVob.*evConst_nu_pcVob)==1;

feas_intersSpEvk = (feas_Evk.*feas_Spont)==1; %the intersection of both feasible sets
feas_unnSpkEvk=(feas_Evk+feas_Spont)>=1; %the union of both feasible sets


pctEv_R_pcVob=sum(feas_reg)/tot_gs
pctEv_Cov_pcVob=sum(evConst_Cov_pcVob)/tot_gs
pctEv_corr_pcVob=sum(evConst_corr_pcVob)/tot_gs
pctEv_Vr=sum(evConst_Vr_pcVob)/tot_gs
pctEv_all=sum(feas_Evk)/tot_gs

pctBoth_inter=sum(feas_intersSpEvk)/tot_gs
pause


%--- Constraint within PC: nu_ev>nu_sp  ------
bConst_nu_pc = mn_pc_frEvk > mn_pc_frSpon;
%--- Constraint within PC: rho_ev<rho_sp  ------
bConst_rho_pc = mnCorr_pc_frEvk < mnCorr_pc_frSpon ;
%--- Constraint within PC: fano_ev<fano_sp  ------
bConst_fano_pc = mnFano_pc_Evk<mnFano_pc_Spon;

%--- Constraint within OB: nu_ev>nu_sp  ------
bConst_nu_ob = mn_ob_frEvk > mn_ob_frSpon;
%--- Constraint within OB: Var_ev>Var_sp  ------
bConst_Vr_ob = vr_ob_frEvk>vr_ob_frSpon;

feas_stCh = (bConst_nu_pc.*bConst_rho_pc.*bConst_fano_pc.*bConst_nu_ob.*bConst_Vr_ob)==1;

finPrm=feas_intersSpEvk.*feas_stCh;

pctNuPC=sum(bConst_nu_pc)/tot_gs
pctNuOB=sum(bConst_nu_ob)/tot_gs
pctRhoPC=sum(bConst_rho_pc)/tot_gs
pctFanoPC=sum(bConst_fano_pc)/tot_gs
pctVrOB=sum(bConst_Vr_ob)/tot_gs
pctStCh=sum(feas_stCh)/tot_gs

pct_all=sum(finPrm)/tot_gs


save dBothSpEv_st finPrm gE* gI*

%bar plots to show percentage of parms satisfy constraints
%se_mean: OB>PC rate and Ev>Spont (combine 4)
%spConst_fano_pcVob: FF of PC>OB in spon; spConst_corr_pcVob:spt PC Corr > OB
%evConst_Cov_pcVob: Ev cov PC<OB; evConst_Vr_pcVob: Ev var PC<OB
%evConst_corr_pcVob: Ev corr PC<OB; bConst_Vr_ob: OB-var-Ev>Sp
%bConst_fano_pc=PC-fano-Ev<Sp; bConst_rho_pc: PC-corr-Ev<Sp
dat_bar_Const=(sum([spConst_nu_pcVob evConst_nu_pcVob bConst_nu_pc bConst_nu_ob spConst_fano_pcVob spConst_corr_pcVob evConst_Cov_pcVob evConst_Vr_pcVob ...
    evConst_corr_pcVob bConst_Vr_ob bConst_fano_pc bConst_rho_pc])./tot_gs)';
dat_bar_Const=[dat_bar_Const 1-dat_bar_Const];
figure
hold on
bar(dat_bar_Const,'stacked')
set(gca,'FontSize',18)
axis off