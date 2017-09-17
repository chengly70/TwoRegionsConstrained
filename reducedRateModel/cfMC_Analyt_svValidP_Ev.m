%Check Monte Carlo of full 6 WC with the analytic theory
%running analytic again BECAUSE want state variables X
%Spontaneous State

load Parms_evoke
load dBothSpEv_st %has finPrm which is the feasible region for both -Spon -Ev_fXMn
crPC=0.35;


Nc=6;
tau_vec=ones(Nc,1);
rv_vec=s_rv*ones(Nc,1);
sp_vec=s_sp*ones(Nc,1);
CinMat=blkdiag( crOB*ones(3,3)+(1-crOB)*diag(ones(3,1)) , crPC*ones(3,3)+(1-crPC)*diag(ones(3,1)) );
% -- coupling matrix, most entries are fixed --
%Orth: 1(I), 2(M/T), 3(M/T); Retr: 4(I), 5(E), 6(E)
Gm=zeros(Nc,Nc);
g_ieo=0.1;
g_iep=0.1;
%assuming Nc=6
Gm=[0 g_ieo*[1 1] zeros(1,3); zeros(2,Nc); zeros(1,4) g_iep*[1 1]; zeros(2,Nc)];

id1=[];
id2=[];
for j=1:2
    id1=[id1; (1:j)'];
    id2=[id2; (j+1)*ones(j,1)];
end
ind_UpTri_ob=sub2ind([Nc Nc],id1,id2); %indices upper triang Corr/Cov matrix
id1=[];
id2=[];
for j=1:2
    id1=[id1; (1:j)'+3];
    id2=[id2; (j+1)*ones(j,1)+3];
end
ind_UpTri_pc=sub2ind([Nc Nc],id1,id2); %indices upper triang Corr/Cov matrix

len_g=length(gEP); %assuming equal length g's; need len_g^4 total runs
leng_2=len_g^2;
leng_3=len_g^3;

indc_gdPm=find(finPrm);
%rn_i=randperm(length(indc_gdPm));
len_p=length(indc_gdPm);
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

tmSave=100; %how often save the data
for j=1:len_p
    %Dim: (1,gIP) (2,gIO) (3,gEP) (4,gEO)
    [ind_ip, ind_io, ind_ep, ind_eo]=ind2sub([len_g len_g len_g len_g],indc_gdPm(j));
    
    %update coupling matrix
    Gm(1,5:6)=gEP(ind_ep);
    Gm(2:3,1)=gIO(ind_io);
    Gm(4,2:3)=gEO(ind_eo);
    Gm(5:6,4)=gIP(ind_ip);
    
    [convged,Corr_valid,cov_Fa,mn_Fa,cov_Xa,mn_Xa,mean_all]=iter_method(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm,CinMat);
    
    crMat=cov_Fa./sqrt(diag(cov_Fa)*diag(cov_Fa)'); %correl matrix, spiking
    crM_X=cov_Xa./sqrt(diag(cov_Xa)*diag(cov_Xa)'); %correl matrix, activi
    %save analytic
    An_covRt_pc_vec(j,:)=cov_Fa(ind_UpTri_pc);
    An_covRt_ob_vec(j,:)=cov_Fa(ind_UpTri_ob);
    An_vr_F_vec(j,:)=diag(cov_Fa);
    An_mn_F_vec(j,:)=mn_Fa;
    An_CorrRt_ob_vec(j,:)=crMat(ind_UpTri_ob);
    An_CorrRt_pc_vec(j,:)=crMat(ind_UpTri_pc);
    An_covX_pc_vec(j,:)=crM_X(ind_UpTri_pc);
    An_covX_ob_vec(j,:)=crM_X(ind_UpTri_ob);
    An_vr_X_vec(j,:)=diag(cov_Xa);
    An_mn_X_vec(j,:)=mn_Xa;
    
    %monte carlo
    [cov_pcFx,cov_obFx,vr_Fs,mn_Fs,cov_pcX,cov_obX,vr_Xs,mn_Xs]=mc_sixWC(mu_vec,sig_vec,crOB,crPC,nuMax,s_rv,s_sp,Gm);
    
    vrF_matr=sqrt(vr_Fs*vr_Fs'); %6 x 6
    
    
    %save MC
    Mc_covRt_pc_vec(j,:)=cov_pcFx';
    Mc_covRt_ob_vec(j,:)=cov_obFx';
    Mc_vr_F_vec(j,:)=vr_Fs';
    Mc_mn_F_vec(j,:)=mn_Fs';
    Mc_CorrRt_ob_vec(j,:)=(cov_obFx./vrF_matr([7;13;14]))';
    Mc_CorrRt_pc_vec(j,:)=(cov_pcFx./vrF_matr([28;34;35]))';
    Mc_covX_pc_vec(j,:)=cov_pcX';
    Mc_covX_ob_vec(j,:)=cov_obX';
    Mc_vr_X_vec(j,:)=vr_Xs';
    Mc_mn_X_vec(j,:)=mn_Xs';
    
    if(j==tmSave || j==len_p)
        save dcfMcAn_Ev indc_gdPm len_p An_* Mc_* crPC
        tmSave=tmSave+100;
    end

end

