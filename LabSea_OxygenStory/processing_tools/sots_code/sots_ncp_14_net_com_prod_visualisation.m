figure('units','normalized','outerposition',[0 0 1 1])

% Plots cummulative net community production
subplot(3,1,1)
yyaxis left
plot(mooring_data.time,cumsum(mooring_data.ncp_O2_molm2hr./(mooring_data.mld_smooth)),'LineWidth',2)
ylabel('mol of O2 per m^{3}')

yyaxis right
plot(mooring_data.time,mooring_data.ncp_O2_molm2_cumsum,'LineWidth',2)
ylabel('mol of O2 per m^{2}')
title('Net Community Production (cummulative)')
grid on
dim = [.795 .66 .1 .1];
str = strcat('NCP total:',compose("%5.1f",mooring_data.ncp_O2_molm2_cumsum(end)),'mol of O2 m^{-2}');
annotation('textbox',dim,'String',str,'FitBoxToText','on');
xlim([mooring_data.time(1) mooring_data.time(end)]);
datetick('x','mmm','KeepLimits');
set(gca,'Fontsize',12)

% Plots smoothed mixed layer depth
subplot(3,1,2)
ylabel('m')
plot(mooring_data.time,mooring_data.mld_smooth,'LineWidth',2)
title('Mixed Layer Depth (smoothed)')
xlim([mooring_data.time(1) mooring_data.time(end)]);
datetick('x','mmm','KeepLimits');
set(gca,'YDir','reverse')
grid on
set(gca,'Fontsize',12)

% Plots oxygen exchange
subplot(3,1,3)
hold on
set(gca,'Fontsize',12)
plot(mooring_data.time,cumsum(mooring_data.dox2_bubbles_molm2),'LineWidth',2)
plot(mooring_data.time,cumsum(mooring_data.dox2_gas_exchange_molm2),'LineWidth',2)
plot(mooring_data.time,cumsum(mooring_data.dox2_biology_molm2),'LineWidth',2)

if exchange_choice == 1
    
    legend('Bubbles','Air-sea gas exchange','Biology');
    
elseif exchange_choice == 2
    
    plot(mooring_data.time,cumsum(mooring_data.dox2_eddy_diffusion_molm2),'LineWidth',2);
    legend('Bubbles','Air-sea gas exchange','Biology','Eddy diffusion')
    
elseif exchange_choice == 3
    
    plot(mooring_data.time,cumsum(mooring_data.dox2_eddy_diffusion_molm2),'LineWidth',2);
    plot(mooring_data.time,cumsum(mooring_data.dox2_entrainment_molm2),'LineWidth',2);
    legend('Bubbles','Air-sea gas exchange','Biology','Eddy diffusion','Entrainment');
    
elseif exchange_choice == 4
    
    
end

ylabel('Cummulative O2 (mol/m^2)')
title('Measured oxygen calculated contributions')
xlim([mooring_data.time(1) mooring_data.time(end)])
datetick('x','mmm','KeepLimits');
grid on

% Plots measured oxygen, oxygen solubility, and calculated physical oxygen
% if exchange_choice==1
%     
%     set(gca,'Fontsize',12)
%     hold on
%     plot(mooring_data.time,mooring_data.dox2_molm3,'LineWidth',2)
%     plot(mooring_data.time,mooring_data.dox2_sol_molm3,'LineWidth',2)
%     plot(mooring_data.time,mooring_data.dox2_phys_molm3,'LineWidth',2)
%     title('Oxygen record (atmospheric exchange only)')
%     legend('Measured','Solubility','Physical model')
%     ylabel('O2 (mol/m^3)')
%     xlim([mooring_data.time(1) mooring_data.time(end)])
%     datetick('x','mmm','KeepLimits');
%     grid on
% 
% elseif exchange_choice==2
%     
%     set(gca,'Fontsize',12)
%     hold on
%     plot(mooring_data.time,mooring_data.dox2_molm3,'LineWidth',2)
%     plot(mooring_data.time,mooring_data.dox2_sol_molm3,'LineWidth',2)
%     plot(mooring_data.time,mooring_data.dox2_phys_molm3,'LineWidth',2)
%     title('Oxygen record (atmospheric and eddy exchange)')
%     legend('Measured','Solubility','Physical model')
%     ylabel('O2 (mol/m^3)')
%     xlim([mooring_data.time(1) mooring_data.time(end)])
%     datetick('x','mmm','KeepLimits');
%     grid on
%     
%     
% elseif exchange_choice==3
%     
%         set(gca,'Fontsize',12)
%     hold on
%     plot(mooring_data.time,mooring_data.dox2_molm3,'LineWidth',2)
%     plot(mooring_data.time,mooring_data.dox2_sol_molm3,'LineWidth',2)
%     plot(mooring_data.time,mooring_data.dox2_phys_molm3,'LineWidth',2)
%     title('Oxygen record (atmospheric, eddy, and entrainment exchange)')
%     legend('Measured','Solubility','Physical model')
%     ylabel('O2 (mol/m^3)')
%     xlim([mooring_data.time(1) mooring_data.time(end)])
%     datetick('x','mmm','KeepLimits');
%     grid on
%     
% elseif exchange_choice==4
%     
% end



% subplot(2,1,1)
% yyaxis left
% plot(mooring_data.time,cumsum(mooring_data.ncp_O2_umm2hr./(1E6*mooring_data.mld_smooth)),'LineWidth',2)
% ylabel('mol of O2 per m^{3}')
% 
% yyaxis right
% plot(mooring_data.time,mooring_data.ncp_O2_umm2_cumsum./1E6,'LineWidth',2)
% ylabel('mol of O2 per m^{2}')
% title('Net Community Production (cummulative)')
% grid on
% dim = [.795 .66 .1 .1];
% str = strcat('NCP total:',compose("%5.1f",mooring_data.ncp_O2_umm2_cumsum(end)/1E6),'mol of O2 m^{-2}');
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% xlim([mooring_data.time(1) mooring_data.time(end)]);
% datetick('x','mmm','KeepLimits');
% set(gca,'Fontsize',12)
% 
% subplot(2,1,2)
% yyaxis left
% plot(mooring_data.time,mooring_data.ncp_O2_umm2hr./mooring_data.mld_smooth,'LineWidth',1.5)
% ylim([-1000 1000])
% ylabel('umol of O2 per m^{3} per hour')
% 
% yyaxis right
% plot(mooring_data.time,mooring_data.ncp_O2_umm2hr,'LineWidth',1.5)
% ylabel('umol of O2 per m^{2}')
% title('Net Community Production (hourly)')
% ylim([-1E6 1E6])
% grid on
% set(gca,'Fontsize',12)
% xlim([mooring_data.time(1) mooring_data.time(end)]);
% datetick('x','mmm','KeepLimits');
% 
% 


