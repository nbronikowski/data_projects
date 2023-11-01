close all
%% Plot Glider Section 
if ~exist('plots', 'dir')
    mkdir('plots')
end

if strcmp(plot_name,'plot_integrated_chla')
    figure()
    dtime = datetime(timeg,'ConvertFrom','datenum'); 
    hold on

    ax_ylims = [ 7 17];

    % Find transitions into and out of nighttime
    night_starts = find(diff(a_zenith >= 90) == 1) + 1;
    night_ends = find(diff(a_zenith >= 90) == -1) + 1;

    % In case the data starts or ends during nighttime, add the first or last index
    if a_zenith(1) >= 90
        night_starts = [1; night_starts];
    end
    if a_zenith(end) >= 90
        night_ends = [night_ends; numel(dtime)];
    end

    % Shade nighttime periods
    for i = 1:length(night_starts)
        fill([dtime(night_starts(i)), dtime(night_ends(i)), dtime(night_ends(i)), dtime(night_starts(i))], ...
             [ax_ylims(1), ax_ylims(1), ax_ylims(2), ax_ylims(2)], ...
             [0.9, 0.9, 0.9], 'EdgeColor', 'none');
    end

    %% Plot Mixed Layer Depth and Night Times
    h1 = plot(dtime, pmld_dens, '-k');
    ylabel('z_{mld} (\sigma_t - 0.05 kg m^{-3}) (m)');
    xlim([dtime(1) dtime(end)])

    %% Plot Integrated Chlorophyll
    yyaxis('right'); 
    hold on
    h2 = plot(dtime, intg_chla, '-g','LineWidth',2);

    ylabel('integrated Chlorophyll z_{mld} - z_{mld+50} (mg m^{-2})');
    title('Vertically Integrated Chlorophyll and Mixed Layer Depth');
    legend([h1 h2],{'Mixed Layer Depth','Integrated Chlorophyll'},'Location','best');
    grid on;
    xlim([dtime(1) dtime(end)])
    formatplot

    save_figure(gcf, ['./plots/vertically_integrated_chla'], [7.5 4], '.png', '300');
end


if strcmp(plot_name,'event_section')

    idnan=find(~isnan(gdat.temp));
    idnan = idnan(1:50:end);

    
    figure()
    subplot(411); hold on
    imagescn(x,-y,pg_temp);  shading flat
    colormap(gca,cmocean('thermal'));
%     plot(gdat.timeDateNum(idnan),-gdat.pres(idnan),'-','LineWidth',0.5,'Color',[.6 .6 .6])
    xlim([nanmin(timeg) nanmax(timeg)])
    caxis([-1.5 10]);
    ylabel('Depth (m)')
    cb1=colorbar;
    ylabel(cb1,'Temperature (^oC)')
    title('Glider Unit\_334 Data Segment with Chlorophyll Maximum')
    datetick('x','dd-mmm','keeplimits')
    
    
    subplot(412); hold on
    imagescn(x,-y,pg_salt); shading flat
    colormap(gca,cmocean('haline'))
    cb2=colorbar;
%     plot(gdat.timeDateNum(idnan),-gdat.pres(idnan),'-','LineWidth',0.5,'Color',[.6 .6 .6])
    xlim([nanmin(timeg) nanmax(timeg)])
    caxis([31 33.5]);
    ylabel('Depth (m)')
    ylabel(cb2,'Salinity (PSU)')
    datetick('x','dd-mmm','keeplimits')
    
    ax3=subplot(413); hold on
    imagescn(x,-y,pg_chla); shading flat
    colormap(gca,cmocean('algae'))
%     plot(gdat.timeDateNum(idnan),-gdat.pres(idnan),'-','LineWidth',0.5,'Color',[.6 .6 .6])
    xlim([nanmin(timeg) nanmax(timeg)])
%     set(ax3, 'ColorScale', 'log')
    cb3=colorbar;
    caxis(ax3,[0 2]); % Modify according to your data range
%     cb3.Ticks = [0.0001,0.001,0.01 0.1,1]; % Modify according to your data range
%     cb3.TickLabels = {'0.0001','0.001','0.01','0.1','1'};
    
    ylabel('Depth (m)')
    ylabel(cb3,'Chl_a (mg m^{-3})')
    datetick('x','dd-mmm','keeplimits')
    
    subplot(414); hold on
    imagescn(x,-y,pg_bbp700); shading flat
    colormap(gca,cmocean('matter'))
    cb4=colorbar;
%     plot(gdat.timeDateNum(idnan),-gdat.pres(idnan),'-','LineWidth',0.5,'Color',[.6 .6 .6])
    xlim([nanmin(timeg) nanmax(timeg)])
    caxis([0 0.02])
    datetick('x','dd-mmm','keeplimits')
    ylabel('Depth (m)')
    ylabel(cb4,'b_{bp,700nm} (10^{-4} m^{-1})')
    
    set(gcf,'Color','w')
    save_figure(gcf,['./plots/glider_event_section'],[7 5],'.png','300');
end

if strcmp(plot_name,'thomalla_correction')
    
    idnan=find(~isnan(gdat.temp));
    idnan = idnan(1:50:end);

    figure()
    % First subplot with symlog color scale
    subplot(211); hold on
    pcolor(x,-y,pg_chla); shading flat
    colormap(gca,cmocean('algae'))
    cb1=colorbar;
    caxis([0 1])
    plot(timeg,-pmld_dens,'-m')
%     plot(gdat.timeDateNum(idnan),-gdat.pres(idnan),'-','LineWidth',0.5,'Color',[.6 .6 .6])
    xlim([nanmin(timeg) nanmax(timeg)])
    ylabel('Depth (m)')
    ylabel(cb1,'Chlorophyll\_a (mg m^{-3})')
    title('Raw Glider Wetlabs Chlorophyll\_a Concentration','Color','w')
    cb1.TickLength = 0.02;
    cb1.Color = 'w';
    box on
    set(gca,'XColor','w','YColor','w')
    datetick('x','dd-mmm','keeplimits')
    ylim([-50 0])

    % Second subplot with symlog color scale
    ax2=subplot(212); hold on
    pcolor(x,-y,chla_corr); shading flat
    colormap(gca,cmocean('algae'))
    cb2=colorbar;
    caxis([0 1])
    plot(timeg,-pmld_dens,'-m')
%     plot(gdat.timeDateNum(idnan),-gdat.pres(idnan),'-','LineWidth',0.5,'Color',[.6 .6 .6])
    xlim([nanmin(timeg) nanmax(timeg)])
    ylabel('Depth (m)')
    ylabel(cb2,'Chlorophyll\_a (mg m^{-3})')
    datetick('x','dd-mmm','keeplimits')
    title('Chlorophyll\_a Quenching Correction from Thomalla et al. (2017)','Color','w')
    cb2.TickLength = 0.02;
    cb2.Color = 'w';
    box on
    set(gca,'XColor','w','YColor','w')
    set(gcf,'Color',[0 0 0])
    ylim([-50 0])
    save_figure(gcf,['./plots/chlorophyll_correction'],[7 5],'.png','300');
end


if strcmp(plot_name,'plot_CMEMS_chla')    
    load CMEMS_multi_chla.mat
    cm.chla(cm.chla<0)=NaN;
    
    PanFac = 4; % we want 4 columns
    
    ChlaMax_idx = gdat.timeDateNum>=datenum(2023,08,27) & gdat.timeDateNum<=datenum(2023,08,28,23,59,59);
    lon_ChlaMax = nanmean(gdat.lon(ChlaMax_idx));
    lat_ChlaMax = nanmean(gdat.lat(ChlaMax_idx));
    
    t = tiledlayout(PanFac, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for i = 11:length(cm.time)-2
        gl_idx = gdato.timeDateNum>=cm.time(i) & gdato.timeDateNum<=cm.time(i)+1;
        
        % Use nexttile instead of subplot
        ax = nexttile; hold on
        imagescn(cm.lon,cm.lat,log10(cm.chla(:,:,i)'))
        set(gca, 'ColorScale', 'log')
    
        h1=plot(gdato.lon,gdato.lat,'.b');
        h2=plot(gdato.lon(gl_idx),gdato.lat(gl_idx),'.r');
        h3=plot(lon_ChlaMax,lat_ChlaMax,'+','MarkerSize',10,...
            'Color','k','LineWidth',2);
        borders('Canada','facecolor',[.6 .6 .6])
        colormap(ax,cmocean('algae'))
        
        % Log color axis from 0.000001 to 2
        caxis(gca,[0.1 1]); % Modify according to your data range
        title(datestr(cm.time(i),'yyyy-mmm-dd'))
        ylim([min(cm.lat) max(cm.lat)])
        xlim([min(cm.lon) max(cm.lon)])
    
        if i == 16
            leg=legend([h1 h2 h3],{'Full Glider Mission','Matching Glider Track','Location of Maximum'},'location','best');
        end
        
        if i == 22
            cb = colorbar;
            cb.Layout.Tile = 'east'; % To position the colorbar on the east side of the tiled layout
            
            % Define logarithmically spaced tick marks
            cb.Ticks = [0.0001,0.001,0.01 0.1,1]; % Modify according to your data range
            cb.TickLabels = {'0.0001','0.001','0.01','0.1','1'};
            ylabel(cb,'Chlorophyll\_a (mg m^{-3})')
        end
    end
    title(t,'CMEMS Multi Product Ocean Colour and Glider Track');
    set(gcf,'Color','w')
    
    save_figure(gcf,['./plots/CMEMS+glider_track'],[9 12],'.png',300')
end


