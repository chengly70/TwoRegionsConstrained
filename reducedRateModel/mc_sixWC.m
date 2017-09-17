function [cov_pcFx,cov_obFx,vr_Fs,mn_Fs,cov_pcX,cov_obX,vr_Xs,mn_Xs,covPCOB_F,covPCOB_X]=mc_sixWC(mu_vec,sig_vec,crOB,crPC,nuMax,s_rv,s_sp,gMat)
%[cov_pcFx,cov_obFx,vr_Fs,mn_Fs,cov_pcX,cov_obX,vr_Xs,mn_Xs,covPCOB_F,covPCOB_X]=mc_sixWC(mu_vec,sig_vec,crOB,crPC,nuMax,s_rv,s_sp,gMat);
%Simulate full network of 6 WC, coupled/corrNoise, etc; set tau_m=1
%sigmoidal coupling; nuMax/2*(1+tanh((Inp-s_rv)./s_sp))
%gMat is a 6x6 matrix, representing full coupling

rng('shuffle') %seed random number generator

dt=0.01; % in msec
t_end=5000; % in msec
Lt=floor(t_end/dt)+1; %total num time-steps
sq_dt=1/sqrt(dt);

N_relz=300;

%state variables (variables needed to do sim)
xV1=sig_vec(1)*randn(N_relz,1)+mu_vec(1); %# of cells
xV2=sig_vec(2)*randn(N_relz,1)+mu_vec(2); %# of cells
xV3=sig_vec(3)*randn(N_relz,1)+mu_vec(3); %# of cells
xV4=sig_vec(4)*randn(N_relz,1)+mu_vec(4); %# of cells
xV5=sig_vec(5)*randn(N_relz,1)+mu_vec(5); %# of cells
xV6=sig_vec(6)*randn(N_relz,1)+mu_vec(6); %# of cells
xi=zeros(6*N_relz,1); %white noise forcing
xic_ob=0;
xic_pc=0;

%OUTPUTS
mn_Xs=zeros(6,1);  %mean of Xj
vr_Xs=zeros(6,1);  %var of Xj
cov_obX=zeros(3,1); %only 3 OB Covs Xj
cov_pcX=zeros(3,1); %only 3 PC Covs Xj
mn_Fs=zeros(6,1);  %mean of F(Xj)
vr_Fs=zeros(6,1);  %var of F(Xj)
cov_obFx=zeros(3,1); %only 3 OB Covs F(Xj)
cov_pcFx=zeros(3,1); %only 3 PC Covs F(Xj)
covPCOB_F=zeros(3,3); %all 9 cross covs F(Xj)
covPCOB_X=zeros(3,3); %all 9 cross covs Xj

%reset vals
c_o=sqrt(1-crOB);
c_cob=sqrt(crOB);
c_p=sqrt(1-crPC);
c_cpc=sqrt(crPC);
mu_rs1=0;
mu_rs2=0;
mu_rs3=0;
mu_rs4=0;
mu_rs5=0;
mu_rs6=0;
mu_F1=0; %firing rates
mu_F2=0;
mu_F3=0;
mu_F4=0;
mu_F5=0;
mu_F6=0;
%cov parts
covWhole_F=zeros(6,6); %all cross covs F(Xj)
covWhole_X=zeros(6,6); %all cross covs Xj

for j=2:Lt
    xi=randn(6*N_relz,1); %generate N_relz standard normal (uncorr in time)
    xic_ob=randn;
    xic_pc=randn;
    
    xVm=[xV1 xV2 xV3 xV4 xV5 xV6];
    FxMat = nuMax/2*(1+tanh((xVm-s_rv)./s_sp));
    
    %eqns for cells
    xV1=xV1+dt*(-xV1+mu_vec(1)+sq_dt*sig_vec(1)*(c_o*xi(1:N_relz)+c_cob*xic_ob)+FxMat*(gMat(1,:)') );
    xV2=xV2+dt*(-xV2+mu_vec(2)+sq_dt*sig_vec(2)*(c_o*xi(N_relz+1:2*N_relz)+c_cob*xic_ob)+FxMat*(gMat(2,:)') );
    xV3=xV3+dt*(-xV3+mu_vec(3)+sq_dt*sig_vec(3)*(c_o*xi(2*N_relz+1:3*N_relz)+c_cob*xic_ob)+FxMat*(gMat(3,:)') );
    
    xV4=xV4+dt*(-xV4+mu_vec(4)+sq_dt*sig_vec(4)*(c_p*xi(3*N_relz+1:4*N_relz)+c_cpc*xic_pc)+FxMat*(gMat(4,:)') );
    xV5=xV5+dt*(-xV5+mu_vec(5)+sq_dt*sig_vec(5)*(c_p*xi(4*N_relz+1:5*N_relz)+c_cpc*xic_pc)+FxMat*(gMat(5,:)') );
    xV6=xV6+dt*(-xV6+mu_vec(6)+sq_dt*sig_vec(6)*(c_p*xi(5*N_relz+1:6*N_relz)+c_cpc*xic_pc)+FxMat*(gMat(6,:)') );
    
    %keep running sum of stats/dens
    mu_rs1=mu_rs1+xV1;
    mu_rs2=mu_rs2+xV2;
    mu_rs3=mu_rs3+xV3;
    mu_rs4=mu_rs4+xV4;
    mu_rs5=mu_rs5+xV5;
    mu_rs6=mu_rs6+xV6;
    
    mu_F1=mu_F1+FxMat(:,1);
    mu_F2=mu_F2+FxMat(:,2);
    mu_F3=mu_F3+FxMat(:,3);
    mu_F4=mu_F4+FxMat(:,4);
    mu_F5=mu_F5+FxMat(:,5);
    mu_F6=mu_F6+FxMat(:,6);
    
    covWhole_X=(covWhole_X*(j-2) + xVm'*xVm)./(j-1);     %running sum; entire cov matrix
    covWhole_F=(covWhole_F*(j-2) + FxMat'*FxMat)./(j-1); %running sum; entire cov matrix
    
end
%normalize and store stats/dens
covWhole_F=covWhole_F./N_relz;
covWhole_X=covWhole_X./N_relz;
mn_Xs(1)=sum(mu_rs1)./(N_relz*Lt);
mn_Xs(2)=sum(mu_rs2)./(N_relz*Lt);
mn_Xs(3)=sum(mu_rs3)./(N_relz*Lt);
mn_Xs(4)=sum(mu_rs4)./(N_relz*Lt);
mn_Xs(5)=sum(mu_rs5)./(N_relz*Lt);
mn_Xs(6)=sum(mu_rs6)./(N_relz*Lt);

mn_Fs(1)=sum(mu_F1)./(N_relz*Lt);
mn_Fs(2)=sum(mu_F2)./(N_relz*Lt);
mn_Fs(3)=sum(mu_F3)./(N_relz*Lt);
mn_Fs(4)=sum(mu_F4)./(N_relz*Lt);
mn_Fs(5)=sum(mu_F5)./(N_relz*Lt);
mn_Fs(6)=sum(mu_F6)./(N_relz*Lt);

vr_Xs = diag(covWhole_X)-mn_Xs.^2;
vr_Fs = diag(covWhole_F)-mn_Fs.^2;

cov_obX(1)=covWhole_X(1,2)-mn_Xs(1)*mn_Xs(2);
cov_obX(2)=covWhole_X(1,3)-mn_Xs(1)*mn_Xs(3);
cov_obX(3)=covWhole_X(2,3)-mn_Xs(2)*mn_Xs(3);
cov_pcX(1)=covWhole_X(4,5)-mn_Xs(4)*mn_Xs(5);
cov_pcX(2)=covWhole_X(4,6)-mn_Xs(4)*mn_Xs(6);
cov_pcX(3)=covWhole_X(5,6)-mn_Xs(5)*mn_Xs(6);
cov_obFx(1)=covWhole_F(1,2)-mn_Fs(1)*mn_Fs(2);
cov_obFx(2)=covWhole_F(1,3)-mn_Fs(1)*mn_Fs(3);
cov_obFx(3)=covWhole_F(2,3)-mn_Fs(2)*mn_Fs(3);
cov_pcFx(1)=covWhole_F(4,5)-mn_Fs(4)*mn_Fs(5);
cov_pcFx(2)=covWhole_F(4,6)-mn_Fs(4)*mn_Fs(6);
cov_pcFx(3)=covWhole_F(5,6)-mn_Fs(5)*mn_Fs(6);

covPCOB_X=covWhole_X(4:6,1:3)-mn_Xs(4:6)*mn_Xs(1:3)';
covPCOB_F=covWhole_F(4:6,1:3)-mn_Fs(4:6)*mn_Fs(1:3)';