%program to find the set (gEP,gEO,gIP,gIO) values that make avgNu(OB)>avgNu(PC)
%adding anatom connections (E to I in each pop, use g_eio, g_eip FIXED)
%similar to ~/model/find_[].m but NOW saving the 'spike' stats
%Spontaneous state; save in dSp_st

load Parms_evoke
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

gEP=(0.1:.1:2)';
gEO=(0.1:.1:2)';

gIP=(-0.1:-.1:-2)';
gIO=(-0.1:-.1:-2)';

len_g=length(gEP); %assuming equal length g's; need len_g^4 total runs
leng_2=len_g^2;
leng_3=len_g^3;

%--- OUTPUTS ----
gs_valid=zeros(len_g^4,1); %indicates if the set is valid or not (see loop for lookup)
frate_s=zeros(len_g^4,6); %all 6 firing rates
rtVar_s=zeros(len_g^4,6); %all 6 variances of firing rates
obCovRt_s=zeros(len_g^4,3); %3 covar OB firing rates
pcCovRt_s=zeros(len_g^4,3); %3 covar PC firing rates
obCorrRt_s=zeros(len_g^4,3); %3 corr OB firing rates
pcCorrRt_s=zeros(len_g^4,3); %3 corr PC firing rates

gE=zeros(2,1); %way Analy_wc_loop.m is set-up
gI=zeros(2,1); %way Analy_wc_loop.m is set-up

tic
for j=1:len_g^4
    
    %could put ind_[] calc outside of loop IF len_g isn't too large
    ind_eo=mod(ceil(j./leng_3)-1,len_g)+1; %outter index
    ind_ep=mod(ceil(j./leng_2)-1,len_g)+1; %2nd index
    ind_io=mod(ceil(j./len_g)-1,len_g)+1; %2nd-to-last index
    ind_ip=mod(j-1,len_g)+1; %inner-most index
    
    %update coupling matrix
    Gm(1,5:6)=gEP(ind_ep);
    Gm(2:3,1)=gIO(ind_io);
    Gm(4,2:3)=gEO(ind_eo);   
    Gm(5:6,4)=gIP(ind_ip);
    
    [convged,Corr_valid,cov_Fa,mn_Fa,cov_Xa,mn_Xa,mean_all]=iter_method(Nc,mu_vec,sig_vec,tau_vec,rv_vec,sp_vec,Gm,CinMat);
    
    crMat=cov_Fa./sqrt(diag(cov_Fa)*diag(cov_Fa)'); %correl matrix
    
    frate_s(j,:)=mn_Fa;
    rtVar_s(j,:)=diag(cov_Fa);
    obCovRt_s(j,:)=cov_Fa(ind_UpTri_ob);
    pcCovRt_s(j,:)=cov_Fa(ind_UpTri_pc);
    obCorrRt_s(j,:)=crMat(ind_UpTri_ob); %save upper-tri (no diagonal)
    pcCorrRt_s(j,:)=crMat(ind_UpTri_pc); %save upper-tri (no diagonal)
    
    if(convged && Corr_valid && sum(mn_Fa(1:3)) > sum(mn_Fa(4:6)) ) %direct calc
        gs_valid(j)=3;
    elseif(~convged) %doesn't converge; oscillates
        gs_valid(j)=2;
    elseif(~Corr_valid) %corr at limit; throwout
        gs_valid(j)=1;
    end
    
    if(ind_ip==len_g)
        save dEv_st gs_valid gEO gEP gIO gIP frate_s rtVar_s obCovRt_s obCorrRt_s pcCovRt_s pcCorrRt_s
    end
        
 end
toc

