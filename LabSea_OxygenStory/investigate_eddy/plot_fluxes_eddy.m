clear; clc; close all
load('workspace.mat')

saltx = salt_intps;
tempx = temp_intps;
oxyx = oxy_intps; % umol kg-1
densx = dens0; 
mldx = MLD;

[~,id_eddycentre]=max(eta_surf);
sc1 = Xq(1,id_eddycentre);
sc2 = 3+Xq(1,id_eddycentre)*2;

[~,idx1] = min(abs(Xq(1,:)-sc1));
[~,idx2] = min(abs(Xq(1,:)-sc2));

Xq2 = Xq;
Yq2 = Yq;
Xq2(:,[1:idx1,idx2:end])=[];
Yq2(:,[1:idx1,idx2:end])=[];
saltx(:,[1:idx1,idx2:end])=[];
tempx(:,[1:idx1,idx2:end])=[];
oxyx(:,[1:idx1,idx2:end])=[];
densx(:,[1:idx1,idx2:end])=[];
mldx([1:idx1,idx2:end])=[];


Rq = Xq2-Xq2(1,1);

%% HATUN PAPER - BUOYANT EDDY
re = 13.5;  % EDDY

D = 1; % m bins
cp = 4000; % J kg^-1 c^-1 
r_ref = 300;

rho_ref = 1027.73 ;
T_ref = 3.2;
S_ref = 34.84;
O2_ref = 300; % avg umol/kg of labsea

Rq_m = 1000*Rq;
re_m = 1000*re;
r_ref_m = 1000*r_ref;
ire = find(Rq_m(1,:)<=re_m);

for i = 1:length(Yq2(:,1))
     Se(i,1) = 1/pi/(re_m)^2*trapz(Rq_m(1,ire),2*pi*Rq_m(1,ire).*saltx(i,ire)); % 
     Te(i,1) = 1/pi/(re_m)^2*trapz(Rq_m(1,ire),2*pi*Rq_m(1,ire).*tempx(i,ire)); % deg C
     O2e(i,1) = 1/pi/(re_m)^2*trapz(Rq_m(1,ire),2*pi*Rq_m(1,ire).*oxyx(i,ire)); % umol kg-1
    
     Salt_Flux_z(i,1) = (S_ref-nanmean(Se(i,1)))    /S_ref   * (re/r_ref)^2      ; % [] /m
     Heat_Flux_z(i,1) = (nanmean(Te(i,1)) - T_ref)  *rho_ref*cp  * (re/r_ref)^2   ; % J / m3
     O2_Flux_z(i,1)   = (nanmean(O2e(i,1)) - O2_ref)/1000 *rho_ref  * (re/r_ref)^2    ;% mols/m3
 

end


% Uper layer
idz1 = find(Yq2(:,1)==330);
D = 330;
FW_1=(S_ref-nanmean(Se(1:idz1)))  / S_ref * 100 * D * (re/r_ref)^2;
H_1=(nanmean(Te(1:idz1))-T_ref) * cp * rho_ref / 1000 * D * (re/r_ref)^2;
O2_1=(nanmean(O2e(1:idz1))-O2_ref) * rho_ref / 1000 * D * (re/r_ref)^2;




% Middle Layer
idz2 = find(Yq2(:,1)==330);
idz3 = find(Yq2(:,1)==600);
D = 600-330;
FW_2=(S_ref-nanmean(Se(idz1:idz2)))  / S_ref * 100 * D * (re/r_ref)^2;
H_2=(nanmean(Te(idz1:idz2))-T_ref) * cp * rho_ref / 1000 * D * (re/r_ref)^2;
O2_2=(nanmean(O2e(idz1:idz2))-O2_ref) * rho_ref / 1000 * D * (re/r_ref)^2;


% Lower Layer
idz2 = find(Yq2(:,1)==330);
idz3 = find(Yq2(:,1)==1000);
D = 1000-330;
FW_3=(S_ref-nanmean(Se(idz1:idz2)))  / S_ref * 100 * D * (re/r_ref)^2;
H_3=(nanmean(Te(idz1:idz2))-T_ref) * cp * rho_ref / 1000 * D * (re/r_ref)^2;
O2_3=(nanmean(O2e(idz1:idz2))-O2_ref) * rho_ref / 1000 * D * (re/r_ref)^2;


table_Headers = {'FW','Heat','Oxygen'};
table_HeadersUnits = {'cm','kJm-2','molm-2'};
rowNames = {'Surface','Intermediate','Deep'};
vals= [FW_1,H_1,O2_1; ...
       FW_2,H_2,O2_2; ...
       FW_3,H_3,O2_3];
% Combine headers and units for clarity
combinedHeaders = strcat(table_Headers, ' [', table_HeadersUnits, ']');

dattab=array2table(vals, 'VariableNames', combinedHeaders, 'RowNames', rowNames);

% Display the table
disp(dattab);

20*11.36*pi*re_m^2/0.37e12*100

%%

Rq2 = [-fliplr(Rq),Rq(:,2:end)];
Yq22 = [fliplr(Yq2),Yq2(:,2:end)];
salt2 = [fliplr(saltx(:,1:end-1)),(saltx)];
temp2 = [fliplr(tempx(:,1:end-1)),(tempx)];
oxy2  = [fliplr(oxyx(:,1:end-1)),(oxyx)];
dens2  = [fliplr(densx(:,1:end-1)),(densx)];
mld2  = [fliplr(mldx(1:end-1)'),mldx'];

% Eddy Stats
figure()

t=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile; hold on
imagescn(Rq2,-Yq22,salt2);
contour(Rq2,-Yq22,dens2,'m','levellist',[27.68 27.72],'showtext','on');
plot(Rq2(1,:),-mld2,'k','LineWidth',2);
cb1=colorbar; colormap(gca,cmocean('haline',10))
caxis([nanmin(salt2(:))-0.001 nanmax(salt2(:))+0.001])
ylim([-1010, 0])
xlabel('Eddy Radius / km')
ylabel('Depth / m');
ylabel(cb1,'S / PSU')
formatplot
cb1.Location='southoutside';

nexttile; hold on
imagescn(Rq2,-Yq22,temp2);
contour(Rq2,-Yq22,dens2,'m','levellist',[27.68 27.72],'showtext','on');
plot(Rq2(1,:),-mld2,'k','LineWidth',2);
cb2=colorbar; colormap(gca,cmocean('thermal',10))
caxis([nanmin(temp2(:))-0.001 nanmax(temp2(:))+0.001])
ylim([-1010, 0])
xlabel('Eddy Radius / km')
set(gca,'YTickLabel','')
ylabel(cb2,'T / ^oC')
formatplot
cb2.Location='southoutside';

nexttile; hold on
imagescn(Rq2,-Yq22,oxy2)
contour(Rq2,-Yq22,dens2,'m','levellist',[27.68 27.72],'showtext','on');
plot(Rq2(1,:),-mld2,'k','LineWidth',2);
cb3=colorbar; load('odv_cmap.mat');
colormap(gca,cmap(1:20:end,:));
caxis([nanmin(oxy2(:))-0.001 nanmax(oxy2(:))+0.001])
ylim([-1010, 0])
set(gca,'YTickLabel','')
xlabel('Eddy Radius / km')
ylabel(cb3,'O_2 / \mumol kg^{-1}')
formatplot
cb3.Location='southoutside';


nexttile; 
plot(Salt_Flux_z*1000,-Yq2(:,1),'b'); hold on;
plot([0 0],[-1000 0],':k');
xlabel('[mm / m]')
ylabel('Depth / m');
formatplot; 
grid on
hline(-330,':k')
hline(-600,':k')
ylim([-1000, 0])
% title('FW Contribution')

nexttile; 
plot(Heat_Flux_z/1000,-Yq2(:,1),'r'); hold on
% plot([0 0],[-1000 0],':k');
set(gca,'YTickLabel','')
xlabel('[KJ / m^3]');
grid on
formatplot;
hline(-330,':k')
hline(-600,':k')
ylim([-1000, 0])
% title('Heat Contribution')

nexttile;
plot(O2_Flux_z,-Yq2(:,1),'k'); hold on
xlabel('[mol / m^3]');
set(gca,'YTickLabel','')
grid on
hline(-330,':k')
hline(-600,':k')
formatplot
ylim([-1000, 0])
% title('O_2 Contribution')


t.TileSpacing='tight';
subtitle(t,'Pearldiver 7-10 February Eddy: FW, Heat and Oxygen Contributions');

save_figure(gcf,'./../plots/pearldiver_eddy_contribution',[7.5 5],'.png','300')