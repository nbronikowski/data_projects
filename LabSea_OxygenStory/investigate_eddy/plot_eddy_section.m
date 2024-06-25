close all

figure()

t=tiledlayout('flow','TileSpacing','compact','Padding','compact')

nexttile(t,1,[1 3]); hold on
ax1 = gca;
h1=plot(Xq(1,:),eta_surf,'k');
ylim([0.4 0.47])
ylabel('\eta (steric 20m) / m')
[~,id_eddycentre]=max(eta_surf);
h2=vline(Xq(1,id_eddycentre),'k');
h3=vline(3,'--r');
vline(3+Xq(1,id_eddycentre)*2,'--r')
leg=legend([h1 h2 h3],{'Steric Height','Eddy Centre','Eddy Extent'});
title('Glider Track through Eddy (Feb  7-10, 2020)')
xlim([0 55])

nexttile([2 3]); hold on; ax2=gca;
imagescn(Xq,-Yq,oxy_intps); shading flat
plot(x(~isnan(y)),-y(~isnan(y)),':k','linewidth',0.5)
% scatter(x,-y,10,dat.O2mmolkg(id))
load('odv_cmap.mat');
colormap(gca,cmap); cb=colorbar;
ylabel(cb,'O_2 / \mumol kg^{-1}')
ylabel('depth / m'); xlabel('glider track distance / km')
contour(Xq,-Yq,dens0,'levellist',[27.68, 27.72, 27.74],'color','m','showtext','on'); 
plot(Xq(1,:),-MLD,'k','LineWidth',2)
caxis([280 300])
vline(Xq(1,id_eddycentre),'k')
vline(3,'--r')
vline(3+Xq(1,id_eddycentre)*2,'--r')
xlim([0 55])


save_figure(gcf,'./../plots/pearldiver_feb_eddy_section',[7.5 3.5],'.png','300');