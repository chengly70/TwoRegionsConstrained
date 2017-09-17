%script to generate heterog thresholds (ThresOB, ThresPC)
%and connectivity matrices ()
%the results are SAVED in 

Npc=100;
Nob=100;

Npce=Npc*.8; %# excit in PC
Npci=Npc*.2; %# inhib in PC
Nobe=Nob*.2; %# excit in OB; M/T cells
Nobi=Nob*.8; %# inhib in OB; WAY more granule cells

percU_pc=(0.95: -0.9/(Npc-1): 0.05)'; %equally spaced CDF; corresp perct% for lognormal
percU_ob=(0.95: -0.9/(Nob-1): 0.05)'; %equally spaced CDF

sig_th=0.2;

ThresOB=zeros(Nob,1);
ThresPC=zeros(Npc,1);

%----intrinsic heterog-----
ThresOB=logninv(percU_ob,0,sig_th);%exp(sig_th*randn(Nob,1));
ThresPC=logninv(percU_pc,0,sig_th);

%--- connectivity matrices ---
fr_ItoEob=0.3;    %gIO
fr_EtoIob=0.3;    %fractions; connectiv 
fr_ItoIob=0.3;
fr_EtoEob=0.3;     
fr_ItoEpc=0.3;    %gIP
fr_EtoIpc=0.3;    %fractions; connectiv 
fr_ItoIpc=0.3;
fr_EtoEpc=0.3;
fr_obToI=0.3;      %gEO, fract of cross (OB->PC,I)
fr_obToE=0.3;      %fract of cross (OB->PC,E)
fr_pcToI=0.3;      %gEP, fract of cross (PC->OB,I)
fr_pcToE=0.3;     %fract of cross (PC->PC,E)

frRob=0.4; %frac of correl saved; total approx. 
frRpc=0.4;

%___avg number of [] inputs to []___
Nm_eio=round(fr_ItoEob*Nobi);  
Nm_ieo=round(fr_EtoIob*Nobe);  
Nm_iio=round(fr_ItoIob*Nobi);  
Nm_eeo=round(fr_EtoEob*Nobe);  
Nm_eip=round(fr_ItoEpc*Npci);  
Nm_iep=round(fr_EtoIpc*Npce);  
Nm_iip=round(fr_ItoIpc*Npci);  
Nm_eep=round(fr_EtoEpc*Npce);  
Nm_iob=round(fr_obToI*Nobe);
Nm_eob=round(fr_obToE*Nobe);
Nm_ipc=round(fr_pcToI*Npce);
Nm_epc=round(fr_pcToE*Npce);


%deterministic number of connections
wEEo=ones(Nobe,1)*Nm_eeo;%poissrnd(Nm_eeo,Nobe,1);
wEIo=ones(Nobe,1)*Nm_eio;%poissrnd(Nm_eio,Nobe,1);
wIEo=ones(Nobi,1)*Nm_ieo;%poissrnd(Nm_ieo,Nobi,1);
wIIo=ones(Nobi,1)*Nm_iio;%poissrnd(Nm_iio,Nobi,1);
wEEp=ones(Npce,1)*Nm_eep;%poissrnd(Nm_eep,Npce,1);
wEIp=ones(Npce,1)*Nm_eip;%poissrnd(Nm_eip,Npce,1);
wIEp=ones(Npci,1)*Nm_iep;%poissrnd(Nm_iep,Npci,1);
wIIp=ones(Npci,1)*Nm_iip;%poissrnd(Nm_iip,Npci,1);
wiPO=ones(Npci,1)*Nm_iob;%poissrnd(Nm_iob,Npci,1);
wePO=ones(Npce,1)*Nm_eob;%poissrnd(Nm_eob,Npce,1);
wiOP=ones(Nobi,1)*Nm_ipc;%poissrnd(Nm_ipc,Nobi,1);
weOP=ones(Nobe,1)*Nm_epc;%poissrnd(Nm_epc,Nobe,1);


J_oo=zeros(Nob,Nob); %weight matrix
for j=1:Nob
    if(j<=Nobe)
        %ind=randi(Nob,1,wOO(j));  %some repeats
        ind=randperm(Nobe,wEEo(j))';
        J_oo(j,ind)=1;
        ind=randperm(Nobi,wEIo(j))';
        J_oo(j,ind+Nobe)=1; %add E-ind (Nobe)
    else
        ind=randperm(Nobe,wIEo(j-Nobe))';
        J_oo(j,ind)=1;
        ind=randperm(Nobi,wIIo(j-Nobe))';
        J_oo(j,ind+Nobe)=1; %add E-ind (Nobe)
    end
end
%J_oo=ones(Nob,Nob); %all-to-all, scale 4 blocks by each Nm
J_oo=J_oo-diag(diag(J_oo)); %remove autaptic coupling

J_pp=zeros(Npc,Npc);
for j=1:Npc
    if(j<=Npce)
        %ind=randi(Nob,1,wOO(j));  %some repeats
        ind=randperm(Npce,wEEp(j))';
        J_pp(j,ind)=1;
        ind=randperm(Npci,wEIp(j))';
        J_pp(j,ind+Npce)=1; %add E-ind (Npce)
    else
        ind=randperm(Npce,wIEp(j-Npce))';
        J_pp(j,ind)=1;
        ind=randperm(Npci,wIIp(j-Npce))';
        J_pp(j,ind+Npce)=1; %add E-ind (Npce)
    end
end
%J_pp=ones(Npc,Npc); %all-to-all, scale 4 blocks by each Nm
J_pp=J_pp-diag(diag(J_pp)); %remove autaptic coupling

J_op=zeros(Nob,Npc);
for j=1:Nob
   if(j<=Nobe)
       ind=randperm(Npce,weOP(j))';
       J_op(j,ind)=1;
   else
       ind=randperm(Npce,wiOP(j-Nobe))';
       J_op(j,ind)=1;
   end
end
%J_op=ones(Nob,Npc); %all-to-all

J_po=zeros(Npc,Nob); %weight matrix
for j=1:Npc
   if(j<=Npce)
      ind=randperm(Nobe,wePO(j))';
      J_po(j,ind)=1;
   else
       ind=randperm(Nobe,wiPO(j-Npce))';
       J_po(j,ind)=1;
   end
end
%J_po=ones(Npc,Nob); %all-to-all

%--trimming so don't have such large connectivity martrices--
nm_cols=max(sum(J_oo,2));
W_oo=zeros(Nob,nm_cols);
nm_cols=max(sum(J_pp,2));
W_pp=zeros(Npc,nm_cols);
nm_cols=max(sum(J_op,2));
W_op=zeros(Nob,nm_cols);
nm_cols=max(sum(J_po,2));
W_po=zeros(Npc,nm_cols);
for j=1:Npc
   indx=find(J_pp(j,:));
   W_pp(j,1:length(indx))=indx;
   indx=find(J_po(j,:));
   W_po(j,1:length(indx))=indx;
end
for j=1:Nob
   indx=find(J_oo(j,:));
   W_oo(j,1:length(indx))=indx;
   indx=find(J_op(j,:));
   W_op(j,1:length(indx))=indx;
end
Nconn_oo=W_oo; %store Nconn_ei to keep track of actual number of connections
Nconn_pp=W_pp;
Nconn_op=W_op;
Nconn_po=W_po;
Nconn_oo(Nconn_oo~=0)=1;
Nconn_pp(Nconn_pp~=0)=1;
Nconn_op(Nconn_op~=0)=1;
Nconn_po(Nconn_po~=0)=1;
%in C: 0,1,..,N-1 and -1 is where stop in C for-loop
W_oo=W_oo-1;
W_pp=W_pp-1;
W_op=W_op-1;
W_po=W_po-1;

%choose random rho_OB pairs
mt_ee=[];
for j=1:(Nob-1)
    mt_ee=[mt_ee; (j-1)*ones(Nob-j,1) (j:Nob-1)'];
end
r_inde=unique( sort(randi(Nob*(Nob-1)/2,round(frRob*Nob*(Nob-1)/2),1)) );
szRob=length(r_inde); %total number of rho_OB pairs
id1_ob=mt_ee(r_inde,1);
id2_ob=mt_ee(r_inde,2);

%choose random rho_PC pairs
mt_ee=[];
for j=1:(Npc-1)
    mt_ee=[mt_ee; (j-1)*ones(Npc-j,1) (j:Npc-1)'];
end
r_inde=unique( sort(randi(Npc*(Npc-1)/2,round(frRpc*Npc*(Npc-1)/2),1)) );
szRpc=length(r_inde); %total number of rho_PC pairs
id1_pc=mt_ee(r_inde,1);
id2_pc=mt_ee(r_inde,2);


save netParms ThresOB ThresPC sig_th percU_* Npc Nob W_* J_* Nconn_* id1_* id2_* szRpc szRob fr_*
