%displays the firing rate for ALL (combined) data

load dFrate_all

pc_bin=(0: 8/30 : ceil(max(nuPC_ev)))';
figure
hold on
hist(nuPC_sp,pc_bin)
histf(nuPC_ev,pc_bin,'Facecolor','none','Edgecolor','r','LineWidth',3)
set(gca,'FontSize',18)
xlabel('Firing Rate (Hz)')
legend('Spontaneous','Evoked')
title('Piriform Cortex')
dim = [.5 .3 .3 .3];
str1 = ['Spont. Rate: ',num2str(mean(nuPC_sp),2),'\pm',num2str(std(nuPC_sp),2)];
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
dim = [.5 .2 .3 .3];
str2 = ['Evoked Rate: ',num2str(mean(nuPC_ev),2),'\pm',num2str(std(nuPC_ev),2)];
tb=annotation('textbox',dim,'String',str2,'FitBoxToText','on');
set(tb,'FontSize',18)

ob_bin=(0: 1 : ceil(max(nuOB_ev))-1)' + 0.5;
figure
hold on
hist(nuOB_sp,ob_bin)
histf(nuOB_ev,ob_bin,'Facecolor','none','Edgecolor','r','LineWidth',3)
set(gca,'FontSize',18)
xlabel('Firing Rate (Hz)')
legend('Spontaneous','Evoked')
title('Olfactory Bulb')
dim = [.5 .3 .3 .3];
str1 = ['Spont. Rate: ',num2str(mean(nuOB_sp),2),'\pm',num2str(std(nuOB_sp),2)];
tb=annotation('textbox',dim,'String',str1,'FitBoxToText','on');
set(tb,'FontSize',18)
dim = [.5 .2 .3 .3];
str2 = ['Evoked Rate: ',num2str(mean(nuOB_ev),2),'\pm',num2str(std(nuOB_ev),2)];
tb=annotation('textbox',dim,'String',str2,'FitBoxToText','on');
set(tb,'FontSize',18)
