%script to check that the Monte Carlo sims also satisfy the constraints
%contains plots to show g valus & percentage each constraint narrows to,
%set the value in tot_gs below 

load dBothSpEv_st
tot_gs=20^4; %!!!! only count amongst 'valid' models where method converges !!!!
len_g=length(gEO); %assume same length

load dcfMcAn_Sp

%Spontaneous state: calc the mean statistics (across 3) from Monte Carlo sims
mc_pc_frSpon=mean(Mc_mn_F_vec(:,4:6),2);
mc_ob_frSpon=mean(Mc_mn_F_vec(:,1:3),2);
mc_pc_fanoSpon=mean(Mc_vr_F_vec(:,4:6)./Mc_mn_F_vec(:,4:6),2);
mc_ob_fanoSpon=mean(Mc_vr_F_vec(:,1:3)./Mc_mn_F_vec(:,1:3),2);
mc_pc_VrSpon=mean(Mc_vr_F_vec(:,4:6),2);
mc_ob_VrSpon=mean(Mc_vr_F_vec(:,1:3),2);
mcCorr_pc_Spon=mean(Mc_CorrRt_pc_vec,2);
mcCorr_ob_Spon=mean(Mc_CorrRt_ob_vec,2);

%--- Constraint: nu PC < nu OB in Spont  ------
constrMCsp_nu_pcVob = (mc_pc_frSpon < mc_ob_frSpon)==1;

%--- Constraint: FF of PC > OB in Spont  ------
constrMCsp_fano_pcVob = (mc_pc_fanoSpon > mc_ob_fanoSpon)==1; 

%--- Constraint: Correl of PC > OB in Spont  ------
constrMCsp_corr_pcVob = (mcCorr_pc_Spon > mcCorr_ob_Spon)==1;


% now for the Evoked state
load dcfMcAn_Ev


%repeat for Evoked: calc the mean statistics (across 3) from Monte Carlo sims
mc_pc_frEvk=mean(Mc_mn_F_vec(:,4:6),2);
mc_ob_frEvk=mean(Mc_mn_F_vec(:,1:3),2);
mc_pc_VrEvk=mean(Mc_vr_F_vec(:,4:6),2);
mc_ob_VrEvk=mean(Mc_vr_F_vec(:,1:3),2);
mcCov_pc_Evk=mean(Mc_covRt_pc_vec,2);
mcCov_ob_Evk=mean(Mc_covRt_ob_vec,2);
mcCorr_pc_Evk=mean(Mc_CorrRt_pc_vec,2);
mcCorr_ob_Evk=mean(Mc_CorrRt_ob_vec,2);
mc_pc_fanoEvk=mean(Mc_vr_F_vec(:,4:6)./Mc_mn_F_vec(:,4:6),2);%needed for checking state-change

%--- Constraint: nu PC < nu OB in Evoked  ------
constrMCev_nu_pcVob = (mc_pc_frEvk < mc_ob_frEvk)==1;

%--- Constraint: Cov of PC < OB in Evoked  ------
constrMCev_Cov_pcVob = (mcCov_pc_Evk < mcCov_ob_Evk)==1; 

%--- Constraint: Var of PC < OB in Evoked  ------
constrMCev_Vr_pcVob = (mc_pc_VrEvk < mc_ob_VrEvk)==1; 

%--- Constraint: Correl of PC < OB in Evoked  ------
constrMCev_corr_pcVob = (mcCorr_pc_Evk < mcCorr_ob_Evk)==1;

% now for change in BOTH states
%--- Constraint within OB: nu_ev>nu_sp  ------
BconsMC_nu_ob = (mc_ob_frEvk > mc_ob_frSpon)==1;
%--- Constraint within PC: nu_ev>nu_sp  ------
BconsMC_nu_pc = (mc_pc_frEvk > mc_pc_frSpon)==1;
%--- Constraint within PC: rho_ev<rho_sp  ------
BconsMC_rho_pc = (mcCorr_pc_Evk < mcCorr_pc_Spon)==1;
%--- Constraint within OB: Var_ev>Var_sp  ------
BconsMC_Vr_ob = (mc_ob_VrEvk > mc_ob_VrSpon)==1;
%--- Constraint within PC: fano_ev<fano_sp  ------
BconsMC_fano_pc = (mc_pc_fanoEvk < mc_pc_fanoSpon)==1;

%collapse all 4 firing rate constraints to 1
se_mean=(constrMCev_nu_pcVob.*constrMCsp_nu_pcVob.*BconsMC_nu_ob.*BconsMC_nu_pc)==1;

Final_Parms=constrMCsp_nu_pcVob.*constrMCsp_fano_pcVob.*constrMCsp_corr_pcVob.*constrMCev_nu_pcVob...
    .*constrMCev_Cov_pcVob.*constrMCev_Vr_pcVob.*constrMCev_corr_pcVob.*BconsMC_nu_ob.*BconsMC_nu_pc...
    .*BconsMC_rho_pc.*BconsMC_Vr_ob.*BconsMC_fano_pc;


indc_An=find(finPrm); %if only want to use WC reduct

% IF you want to look at those parameter sets that satisfy Monte Carlo as well
indc_gdPm=find(Final_Parms);

indc_badParm=ones(len_g^4,1);

% IF only analytic
len_p = length(indc_An);


gValid_ip=zeros(len_p,1);
gValid_io=zeros(len_p,1);
gValid_ep=zeros(len_p,1);
gValid_eo=zeros(len_p,1);

for j=1:len_p
    
    %Dim: (1,gIP) (2,gIO) (3,gEP) (4,gEO)
    
    % ONLY analytic
    indc_badParm(indc_An(j))=0;
    [ind_ip, ind_io, ind_ep, ind_eo]=ind2sub([len_g len_g len_g len_g],indc_An(j));
    
    %Analytic AND MC
    %indc_badParm(indc_An(indc_gdPm(j)))=0;
    %[ind_ip, ind_io, ind_ep, ind_eo]=ind2sub([len_g len_g len_g len_g],indc_An(indc_gdPm(j)));
    
    gValid_ip(j)=gIP(ind_ip);
    gValid_io(j)=gIO(ind_io);
    gValid_ep(j)=gEP(ind_ep);
    gValid_eo(j)=gEO(ind_eo);
end



% -- check structure for ---
Gv_mat=[gValid_ip gValid_io gValid_ep gValid_eo];

[Um,Sgv,Vm]=svd(Gv_mat,0);

%SVD: subtract mean off first
meanGv = mean(Gv_mat);
[Um2,Sgv2,Vm2]=svd(Gv_mat-repmat(meanGv,length(Gv_mat),1),0);
S_re2=Sgv2;
S_re2(2,2)=0; 
S_re2(3,3)=0;
S_re2(4,4)=0;
Gred_an2=Um2*S_re2*Vm2' + repmat(meanGv,length(Gv_mat),1);


disp(['mean gIP = ',num2str(meanGv(1))])
disp(['mean gIO = ',num2str(meanGv(2))])
disp(['mean gEP = ',num2str(meanGv(3))])
disp(['mean gEO = ',num2str(meanGv(4))])

disp(['EigVec Largest SV = ',num2str(Vm2(:,1)')])

disp(['EigVec 2nd Largest SV = ',num2str(Vm2(:,2)')])

disp(['Perc of var SVD = ',num2str((Sgv(1)+Sgv(2,2))/sum(diag(Sgv)))])

save dBothSt_MC Final_Parms

%bar plots to show percentage of parms satisfy constraints
%se_mean: OB>PC rate and Ev>Spont (combine 4)
%constrMCsp_fano_pcVob: FF of PC>OB in spon; constrMCsp_corr_pcVob:spt PC Corr > OB
%constrMCev_Cov_pcVob: Ev cov PC<OB; constrMCev_Vr_pcVob: Ev var PC<OB
%constrMCev_corr_pcVob: Ev corr PC<OB; BconsMC_Vr_ob: OB-var-Ev>Sp
%BconsMC_fano_pc=PC-fano-Ev<Sp; BconsMC_rho_pc: PC-corr-Ev<Sp
% dat_bar_Const=(sum([se_mean constrMCsp_fano_pcVob constrMCsp_corr_pcVob constrMCev_Cov_pcVob constrMCev_Vr_pcVob ...
%     constrMCev_corr_pcVob BconsMC_Vr_ob BconsMC_fano_pc BconsMC_rho_pc])./tot_gs)';
dat_bar_Const=(sum([constrMCsp_nu_pcVob constrMCev_nu_pcVob BconsMC_nu_pc BconsMC_nu_ob constrMCsp_fano_pcVob constrMCsp_corr_pcVob constrMCev_Cov_pcVob constrMCev_Vr_pcVob ...
     constrMCev_corr_pcVob BconsMC_Vr_ob BconsMC_fano_pc BconsMC_rho_pc])./tot_gs)';
dat_bar_Const=[dat_bar_Const 1-dat_bar_Const];
figure
hold on
bar(dat_bar_Const,'stacked')
set(gca,'FontSize',18)
set(gca,'YLim',[0 .012])
axis off