clear; close all; clc;
clear; clc;
load ADCP_inversion_result.mat

%% Glider DAC and ADCP DAC
figure()
idnan = ~isnan(out.ad2cp_u_dac);
[pf_u,pfs_u]=polyfit(out.ad2cp_u_dac(idnan),out.glider_u_dac(idnan),1);

u_R2 = round(1 - (pfs_u.normr/norm(out.glider_u_dac(idnan) - mean(out.glider_u_dac(idnan))))^2,2)
u_R2_1to1 = round(1 - sum((out.glider_u_dac(idnan) - out.ad2cp_u_dac(idnan)*1).^2)...
    /sum((out.glider_u_dac(idnan) - mean(out.glider_u_dac(idnan))).^2),2);


[pf_v,pfs_v]=polyfit(out.ad2cp_v_dac(idnan),out.glider_v_dac(idnan),1);
v_R2 = round(1 - (pfs_v.normr/norm(out.glider_v_dac(idnan) - mean(out.glider_v_dac(idnan))))^2,2);
v_R2_1to1 = round(1 - sum((out.glider_v_dac(idnan) - out.ad2cp_v_dac(idnan)*1).^2)...
    /sum((out.glider_v_dac(idnan) - mean(out.glider_v_dac(idnan))).^2),2);


subplot(121); hold on
h1=plot(out.ad2cp_u_dac,out.glider_u_dac,'Marker','square',...
    'LineStyle','none','MarkerEdgeColor','k');
h2=plot([-0.5 0.5],[-0.5 0.5],'-','Color',[.6 .6 .6]);
h3=plot([-0.5 0.5],polyval(pf_u,[-0.5 0.5]),'LineStyle',':','Color','r');
legend([h2 h3],{['1:1 line, R^2 = ',num2str(u_R2_1to1)],['Best Fit, R^2 = ',num2str(u_R2)]},...
    'Location','NW')
set(gca,'XGrid','on','YGrid','on')
ylabel('u_{g} (m s^{-1})')
xlabel('u_{adcp} (m s^{-1})')
xlim([-0.35 0.25])
ylim([-0.35 0.25])
formatplot

subplot(122); hold on
h4=plot(out.ad2cp_v_dac,out.glider_v_dac,'Marker','square',...
    'LineStyle','none','MarkerEdgeColor','k');
h5=plot([-0.5 0.5],[-0.5 0.5],'-','Color',[.6 .6 .6]);
h6=plot([-0.5 0.5],polyval(pf_v,[-0.5 0.5]),'LineStyle',':','Color','r');
set(gca,'XGrid','on','YGrid','on')
legend([h5 h6],{['1:1 line, R^2 = ',num2str(v_R2_1to1)],['Best Fit, R^2 = ',num2str(v_R2)]},...
    'Location','NW')
ylabel('v_{g} (m s^{-1})')
xlabel('v_{adcp} (m s^{-1})')
xlim([-0.35 0.25])
ylim([-0.35 0.25])
formatplot

save_figure(gcf,'./plots/ADCP_glider_comparison',[7.5 3.5],'.png','300')


