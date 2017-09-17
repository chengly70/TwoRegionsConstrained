
%plotting comparison of MC and Analytic Theory WITH errobars!

load dcfMcAnEB_rndSmpl_Sp %100 random samples
p_EB=0.95; %Percent confidence in error bars
sc_std=norminv((1+p_EB)/2, 0, 1); %assuming 0<= p_EB < 1

figure
hold on
for j=1:3
    plot(An_CorrRt_pc_vec(:,j),Mc_CorrRt_pc_vec(:,j),'*','color',[0 0.5 0 ],'MarkerSize',12)
    plot(An_CorrRt_ob_vec(:,j),Mc_CorrRt_ob_vec(:,j),'.','MarkerSize',24)
end
mn_diag=min(min([An_CorrRt_ob_vec An_CorrRt_pc_vec Mc_CorrRt_pc_vec Mc_CorrRt_ob_vec]));
mx_diag=max(max([An_CorrRt_ob_vec An_CorrRt_pc_vec Mc_CorrRt_pc_vec Mc_CorrRt_ob_vec]));
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k') %diagonal line
set(gca,'FontSize',18)
xlabel('Analytic')
ylabel('Monte Carlo')
title('Correlation of Firing Rate')

figure
hold on
for j=1:3
    errorbar(An_covRt_pc_vec(:,j),Mc_covRt_pc_vec(:,j),sc_std*sqrt(abs(EBMc_covRt_pc_vec(:,j))),'*','color',[0 0.5 0 ],'MarkerSize',12)
    errorbar(An_covRt_ob_vec(:,j),Mc_covRt_ob_vec(:,j),sc_std*sqrt(abs(EBMc_covRt_ob_vec(:,j))),'.','MarkerSize',24)
end
mn_diag=min(min([An_covRt_ob_vec An_covRt_pc_vec Mc_covRt_pc_vec Mc_covRt_ob_vec]));
mx_diag=max(max([An_covRt_ob_vec An_covRt_pc_vec Mc_covRt_pc_vec Mc_covRt_ob_vec]));
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k') %diagonal line
set(gca,'FontSize',18)
xlabel('Analytic')
ylabel('Monte Carlo')
title('Covariance of Firing Rate')


figure
hold on
for j=1:6
    if(j<=3)
        errorbar(An_vr_F_vec(:,j),Mc_vr_F_vec(:,j),sc_std*sqrt(abs(EBMc_vr_F_vec(:,j))),'.','MarkerSize',24)
    else
        errorbar(An_vr_F_vec(:,j),Mc_vr_F_vec(:,j),sc_std*sqrt(abs(EBMc_vr_F_vec(:,j))),'*','color',[0 0.5 0 ],'MarkerSize',12)
    end
end
mn_diag=min(min([An_vr_F_vec An_vr_F_vec Mc_vr_F_vec Mc_vr_F_vec]));
mx_diag=max(max([An_vr_F_vec An_vr_F_vec Mc_vr_F_vec Mc_vr_F_vec]));
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k') %diagonal line
set(gca,'FontSize',18)
xlabel('Analytic')
ylabel('Monte Carlo')
title('Variance of Firing Rate')

figure
hold on
for j=1:6
    if(j<=3)
        errorbar(An_mn_F_vec(:,j),Mc_mn_F_vec(:,j),sc_std*sqrt(abs(EBMc_mn_F_vec(:,j))),'.','MarkerSize',24)
    else
        errorbar(An_mn_F_vec(:,j),Mc_mn_F_vec(:,j),sc_std*sqrt(abs(EBMc_mn_F_vec(:,j))),'*','color',[0 0.5 0 ],'MarkerSize',12)
    end
end
mn_diag=min(min([An_mn_F_vec An_mn_F_vec Mc_mn_F_vec Mc_mn_F_vec]));
mx_diag=max(max([An_mn_F_vec An_mn_F_vec Mc_mn_F_vec Mc_mn_F_vec]));
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k') %diagonal line
set(gca,'FontSize',18)
xlabel('Analytic')
ylabel('Monte Carlo')
title('Mean Firing Rate')


%--- below is same as above BUT for X (not F), activity (instead of Firing Rate)

figure
hold on
for j=1:3
    errorbar(An_covX_pc_vec(:,j),Mc_covX_pc_vec(:,j),sc_std*sqrt(abs(EBMc_covX_pc_vec(:,j))),'*','color',[0 0.5 0 ],'MarkerSize',12)
    errorbar(An_covX_ob_vec(:,j),Mc_covX_ob_vec(:,j),sc_std*sqrt(abs(EBMc_covX_ob_vec(:,j))),'.','MarkerSize',24)
end
mn_diag=min(min([An_covX_ob_vec An_covX_pc_vec Mc_covX_pc_vec Mc_covX_ob_vec]));
mx_diag=max(max([An_covX_ob_vec An_covX_pc_vec Mc_covX_pc_vec Mc_covX_ob_vec]));
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k') %diagonal line
set(gca,'FontSize',18)
xlabel('Analytic')
ylabel('Monte Carlo')
title('Covariance of X')


figure
hold on
for j=1:6
    if(j<=3)
        errorbar(An_vr_X_vec(:,j),Mc_vr_X_vec(:,j),sc_std*sqrt(abs(EBMc_vr_X_vec(:,j))),'.','MarkerSize',24)
    else
        errorbar(An_vr_X_vec(:,j),Mc_vr_X_vec(:,j),sc_std*sqrt(abs(EBMc_vr_X_vec(:,j))),'*','color',[0 0.5 0 ],'MarkerSize',12)
    end
end
mn_diag=min(min([An_vr_X_vec An_vr_X_vec Mc_vr_X_vec Mc_vr_X_vec]));
mx_diag=max(max([An_vr_X_vec An_vr_X_vec Mc_vr_X_vec Mc_vr_X_vec]));
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k') %diagonal line
set(gca,'FontSize',18)
xlabel('Analytic')
ylabel('Monte Carlo')
title('Variance of X')

figure
hold on
for j=1:6
    if(j<=3)
        errorbar(An_mn_X_vec(:,j),Mc_mn_X_vec(:,j),sc_std*sqrt(abs(EBMc_mn_X_vec(:,j))),'.','MarkerSize',24)
    else
        errorbar(An_mn_X_vec(:,j),Mc_mn_X_vec(:,j),sc_std*sqrt(abs(EBMc_mn_X_vec(:,j))),'*','color',[0 0.5 0 ],'MarkerSize',12)
    end
end
mn_diag=min(min([An_mn_X_vec An_mn_X_vec Mc_mn_X_vec Mc_mn_X_vec]));
mx_diag=max(max([An_mn_X_vec An_mn_X_vec Mc_mn_X_vec Mc_mn_X_vec]));
xl=(mn_diag: (mx_diag-mn_diag)/99: mx_diag)';
plot(xl,xl,'k') %diagonal line
set(gca,'FontSize',18)
xlabel('Analytic')
ylabel('Monte Carlo')
title('Mean of X')
