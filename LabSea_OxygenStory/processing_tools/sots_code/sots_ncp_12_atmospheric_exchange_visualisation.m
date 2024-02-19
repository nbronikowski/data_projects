% Plots the atmospheric gas exchange data
    
figure('units','normalized','outerposition',[0 0 1 1])

% Plots calculated nitrogen and nitrogen solubility
subplot(2,2,1)
set(gca,'Fontsize',12)
hold on
plot(mooring_data.time,mooring_data.N2_molm3)
ylabel('N2 (mol/m^3)')
plot(mooring_data.time,(mooring_data.atmosphericpress_Pa/constants.atm_in_Pa).*mooring_data.N2sol_molm3)
legend('Measured','Solubility')
title('Nitrogen record')
xlim([mooring_data.time(1) mooring_data.time(end)])
datetick('x','mmm','KeepLimits');
grid on

% Plots nitrogen exchange
subplot(2,2,3)
hold on
set(gca,'Fontsize',12)
plot(mooring_data.time,cumsum(mooring_data.N2_bubbles_molm3))
plot(mooring_data.time,cumsum(mooring_data.N2_gas_ex_molm3))
legend('Bubbles','Air-sea gas exchange')
ylabel('Cummulative N2 (mol/m^3)')
title('Measured Nitrogen record components')
xlim([mooring_data.time(1) mooring_data.time(end)])
datetick('x','mmm','KeepLimits');
grid on

% Plots measured oxygen, solubility, and calculated physical oxygen
subplot(2,2,2)
set(gca,'Fontsize',12)
hold on
plot(mooring_data.time,mooring_data.dox2_molm3)
plot(mooring_data.time,(mooring_data.atmosphericpress_Pa/constants.atm_in_Pa).*mooring_data.dox2_sol_molm3)
plot(mooring_data.time,mooring_data.dox2_phys_molm3)
plot(mooring_data.time,mooring_data.sub_mld_dox2_molm3)
legend('Measured','Solubility','Physical model','Sub MLD')
ylabel('O2 (mol/m^3)')
title('Oxygen record (Atmospheric, eddy and entrainment exchange)')
xlim([mooring_data.time(1) mooring_data.time(end)])
datetick('x','mmm','KeepLimits');
grid on

% Plots oxygen exchange
subplot(2,2,4)
hold on
set(gca,'Fontsize',12)
plot(mooring_data.time,cumsum(mooring_data.dox2_bubbles_molm2))
plot(mooring_data.time,cumsum(mooring_data.dox2_gas_exchange_molm2))
plot(mooring_data.time,cumsum(mooring_data.dox2_biology_molm2))

if exchange_choice == 1
    
    legend('Bubbles','Air-sea gas exchange','Biology');
    
elseif exchange_choice == 2
    
    plot(mooring_data.time,cumsum(mooring_data.dox2_eddy_diffusion_molm2));
    legend('Bubbles','Air-sea gas exchange','Biology','Eddy diffusion')
    
elseif exchange_choice == 3
    
    plot(mooring_data.time,cumsum(mooring_data.dox2_eddy_diffusion_molm2));
    plot(mooring_data.time,cumsum(mooring_data.dox2_entrainment_molm2));
    legend('Bubbles','Air-sea gas exchange','Biology','Eddy diffusion','Entrainment');
    
elseif exchange_choice == 4
    
    
end

ylabel('Cummulative O2 (mol/m^2)')
title('Measured oxygen calculated contributions')
xlim([mooring_data.time(1) mooring_data.time(end)])
datetick('x','mmm','KeepLimits');
grid on
