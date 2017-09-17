%Check how g's relate to each other


load dBothSpEv_st %has finPrm which is the feasible region for both -Spon -Ev_fXMn

len_g=length(gEP); %assuming equal length g's; need len_g^4 total runs
leng_2=len_g^2;
leng_3=len_g^3;

indc_gdPm=find(finPrm);

len_p=length(indc_gdPm);
gValid_eo=zeros(len_p,1);
gValid_ep=zeros(len_p,1);
gValid_io=zeros(len_p,1);
gValid_ip=zeros(len_p,1);

for j=1:len_p
    %Dim: (1,gIP) (2,gIO) (3,gEP) (4,gEO)
    [ind_ip, ind_io, ind_ep, ind_eo]=ind2sub([len_g len_g len_g len_g],indc_gdPm(j));
    
    gValid_eo(j)=gEO(ind_eo);
    gValid_ep(j)=gEP(ind_ep);
    gValid_io(j)=gIO(ind_io);
    gValid_ip(j)=gIP(ind_ip);

end


close all

figure(1)
hold on
plot(abs(gValid_io),abs(gValid_ip),'b.','MarkerSize',24)

figure(2)
hold on
plot(gValid_eo,gValid_ep,'b.','MarkerSize',24)

figure(3)
hold on
plot(abs(gValid_io),gValid_eo,'b.','MarkerSize',24)

figure(4)
hold on
plot(abs(gValid_ip),gValid_eo,'b.','MarkerSize',24)


figure(5)
hold on
plot(abs(gValid_ip),gValid_ep,'b.','MarkerSize',24)

figure(6)
hold on
plot(abs(gValid_io),gValid_ep,'b.','MarkerSize',24)




% -- now plot the results from MC on top ---

load dBothSt_MC

len_g=length(gEP); %assuming equal length g's; need len_g^4 total runs
leng_2=len_g^2;
leng_3=len_g^3;

indc_gdPm=find(finPrm);
Sub_indc_gdPm=find(Final_Parms); %AFTER run check_Constr_MC, then can re-plot
                                  %NOTE: Final_Parms is of size indc_gdPm==1 (11 or 143)
len_p=length(Sub_indc_gdPm);
gValid_eo=zeros(len_p,1);
gValid_ep=zeros(len_p,1);
gValid_io=zeros(len_p,1);
gValid_ip=zeros(len_p,1);

for j=1:len_p
    %Dim: (1,gIP) (2,gIO) (3,gEP) (4,gEO)
    [ind_ip, ind_io, ind_ep, ind_eo]=ind2sub([len_g len_g len_g len_g],indc_gdPm(Sub_indc_gdPm(j)));
    
    gValid_eo(j)=gEO(ind_eo);
    gValid_ep(j)=gEP(ind_ep);
    gValid_io(j)=gIO(ind_io);
    gValid_ip(j)=gIP(ind_ip);

end


figure(1)
plot(abs(gValid_io),abs(gValid_ip),'r.','MarkerSize',24)
axis([0 2 0 2])
set(gca,'FontSize',18)
xlabel('|gIO|')
ylabel('|gIP|')

figure(2)
plot(gValid_eo,gValid_ep,'r.','MarkerSize',24)
axis([0 2 0 2])
set(gca,'FontSize',18)
xlabel('gEO')
ylabel('gEP')

figure(3)
plot(abs(gValid_io),gValid_eo,'r.','MarkerSize',24)
axis([0 2 0 2])
set(gca,'FontSize',18)
xlabel('|gIO|')
ylabel('gEO')


figure(4)
hold on
plot(abs(gValid_ip),gValid_eo,'r.','MarkerSize',24)
axis([0 2 0 2])
set(gca,'FontSize',18)
xlabel('|gIP|')
ylabel('gEO')


figure(5)
plot(abs(gValid_ip),gValid_ep,'r.','MarkerSize',24)
axis([0 2 0 2])
set(gca,'FontSize',18)
xlabel('|gIP|')
ylabel('gEP')

figure(6)
hold on
plot(abs(gValid_io),gValid_ep,'r.','MarkerSize',24)
axis([0 2 0 2])
set(gca,'FontSize',18)
xlabel('|gIO|')
ylabel('gEP')


