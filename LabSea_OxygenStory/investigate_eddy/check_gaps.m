clear; clc;
path_name = './../mat_files/'; var_name  = 'pearldiver_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = pearldiver; clear pearldiver;

yq = 0:0.5:1020;
[prof_idx, prof_dir] = findProfiles(dat.depth,'inversion',5,'length',100,'stall',50);

length(unique(dat.profile_index))
length(unique(prof_idx))

[T,ux] = pgrid_columns(prof_idx,dat.depth,dat.adjusted_oxygen_concentration,yq);

gl_time = mean_interp(prof_idx,dat.dateNum,ux,'linear');

figure()
imagescn(gl_time,-yq,T);
colormap(cmocean('thermal'))
caxis([280 310])
datetick('x','keeplimits')

t3=datenum(2020,02,08);% eddy 2
t4=datenum(2020,02,10);
id2 = find(dat.dateNum>t3-4 & dat.dateNum<t4+4);

figure()
scatter(dat.dateNum(id2),-dat.depth(id2),20,dat.adjusted_oxygen_concentration(id2),'filled')
colormap(cmocean('thermal'))
caxis([280 310])
datetick('x','keeplimits')

t = dat.dateNum(id2);
t = t-t(1);
[xq,yq]=meshgrid([0:1/24:max(t)],[0:5:1020]);
Zq=barnes(t,dat.depth(id2),dat.adjusted_oxygen_concentration(id2),xq,yq,0.5,20);

figure()
plot(dat.longitude(id2),dat.latitude(id2),'.')


figure()
subplot(211)
scatter(t,-dat.depth(id2),20,dat.adjusted_oxygen_concentration(id2),'filled')
colormap(cmocean('thermal'))
caxis([280 310])

subplot(212)
imagescn(xq,-yq,Zq)
colormap(cmocean('thermal'))
caxis([280 310])