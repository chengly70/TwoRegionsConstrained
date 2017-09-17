%Check Monte Carlo of full 6 WC with the analytic theory; same as 
%cfMC_Analy_rndSmpl but WITH error bars!! (using mx_WCebt.c)

load Parms_spont
load dBothSpEv_pt2 %has finPrm which is the feasible region for both -Spon -Ev_fXMn

len_g=length(gEP); %assuming equal length g's; need len_g^4 total runs
leng_2=len_g^2;
leng_3=len_g^3;

rnd_permL=randperm(len_g^4);
index_sim=rnd_permL(1:100)'; %only do 100 random samples


len_p=length(index_sim);
%-- outputs to save --
An_covRt_pc_vec=zeros(len_p,3);
Mc_covRt_pc_vec=zeros(len_p,3);
An_covRt_ob_vec=zeros(len_p,3);
Mc_covRt_ob_vec=zeros(len_p,3);
An_CorrRt_pc_vec=zeros(len_p,3);
Mc_CorrRt_pc_vec=zeros(len_p,3);
An_CorrRt_ob_vec=zeros(len_p,3);
Mc_CorrRt_ob_vec=zeros(len_p,3);
An_vr_F_vec=zeros(len_p,6);
Mc_vr_F_vec=zeros(len_p,6);
An_mn_F_vec=zeros(len_p,6);
Mc_mn_F_vec=zeros(len_p,6);
An_covX_pc_vec=zeros(len_p,3);
Mc_covX_pc_vec=zeros(len_p,3);
An_covX_ob_vec=zeros(len_p,3);
Mc_covX_ob_vec=zeros(len_p,3);
An_vr_X_vec=zeros(len_p,6);
Mc_vr_X_vec=zeros(len_p,6);
An_mn_X_vec=zeros(len_p,6);
Mc_mn_X_vec=zeros(len_p,6);
EBMc_covRt_pc_vec=zeros(len_p,3);
EBMc_covRt_ob_vec=zeros(len_p,3);
EBMc_vr_F_vec=zeros(len_p,6);
EBMc_mn_F_vec=zeros(len_p,6);
EBMc_covX_pc_vec=zeros(len_p,3);
EBMc_covX_ob_vec=zeros(len_p,3);
EBMc_vr_X_vec=zeros(len_p,6);
EBMc_mn_X_vec=zeros(len_p,6);

tmSave=10; %how often save the data
for j=1:len_p
%Dim: (1,gIP) (2,gIO) (3,gEP) (4,gEO)
[ind_ip, ind_io, ind_ep, ind_eo]=ind2sub([len_g len_g len_g len_g],index_sim(j));
gE=zeros(2,1);
gI=zeros(2,1);
gE(1)=gEO(ind_eo);
gE(2)=gEP(ind_ep);
gI(1)=gIO(ind_io);
gI(2)=gIP(ind_ip);

%analytic calc
[convged,Corr_valid,rate_all,varRate_all,CovRt_ob,CovRt_pc,CorRt_ob,CorRt_pc,mean_all,var_all,cov_ob,cov_pc]=Analy_wc_loop(mu_vec,sig_vec,crOB,crPC,nuMax,s_rv,s_sp,gE,gI,ee_layer); 

gMat=zeros(6,6);
gMat(2:3,1)=gI(1);
gMat(5:6,4)=gI(2);
gMat(1,5:6)=gE(2);
%gMat(2:3,5:6)=ee_layer*gE(2);
gMat(4,2:3)=gE(1);
gMat(5:6,2:3)=ee_layer*gE(1);

%monte carlo, error bars
[cov_pcFx,cov_obFx,vr_Fs,mn_Fs,cov_pcX,cov_obX,vr_Xs,mn_Xs,EBcov_pcFx,EBcov_obFx,EBvr_Fs,EBmn_Fs,EBcov_pcX,EBcov_obX,EBvr_Xs,EBmn_Xs]=mx_WCebt(mu_vec,sig_vec,crOB,crPC,nuMax,s_rv,s_sp,gMat);


vrF_matr=sqrt(vr_Fs*vr_Fs'); %6 x 6

An_covRt_pc_vec(j,:)=CovRt_pc';
Mc_covRt_pc_vec(j,:)=cov_pcFx';
An_covRt_ob_vec(j,:)=CovRt_ob';
Mc_covRt_ob_vec(j,:)=cov_obFx';
An_vr_F_vec(j,:)=varRate_all';
Mc_vr_F_vec(j,:)=vr_Fs';
An_mn_F_vec(j,:)=rate_all';
Mc_mn_F_vec(j,:)=mn_Fs';
An_CorrRt_ob_vec(j,:)=CorRt_ob';
Mc_CorrRt_ob_vec(j,:)=(cov_obFx./vrF_matr([7;13;14]))';
An_CorrRt_pc_vec(j,:)=CorRt_pc';
Mc_CorrRt_pc_vec(j,:)=(cov_pcFx./vrF_matr([28;34;35]))';

An_covX_pc_vec(j,:)=cov_pc(:,end)';
Mc_covX_pc_vec(j,:)=cov_pcX';
An_covX_ob_vec(j,:)=cov_ob(:,end)';
Mc_covX_ob_vec(j,:)=cov_obX';
An_vr_X_vec(j,:)=var_all(:,end)';
Mc_vr_X_vec(j,:)=vr_Xs';
An_mn_X_vec(j,:)=mean_all(:,end)';
Mc_mn_X_vec(j,:)=mn_Xs';

%error bar part
EBMc_covRt_pc_vec(j,:)=EBcov_pcFx';
EBMc_covRt_ob_vec(j,:)=EBcov_obFx';
EBMc_vr_F_vec(j,:)=EBvr_Fs';
EBMc_mn_F_vec(j,:)=EBmn_Fs';
%EBMc_CorrRt_ob_vec(j,:)=(cov_obFx./vrF_matr([7;13;14]))';
%EBMc_CorrRt_pc_vec(j,:)=(cov_pcFx./vrF_matr([28;34;35]))';
EBMc_covX_pc_vec(j,:)=EBcov_pcX';
EBMc_covX_ob_vec(j,:)=EBcov_obX';
EBMc_vr_X_vec(j,:)=EBvr_Xs';
EBMc_mn_X_vec(j,:)=EBmn_Xs';

if(j==tmSave || j==len_p)
    save dcfMcAnEB_rndSmpl_Sp index_sim len_p An_* Mc_* EBMc_*
    tmSave=tmSave+10;
end

end

