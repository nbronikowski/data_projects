clear; clc; close all;

load('azfp_out.mat');
load('glider_ctd_dcast.mat')

%% Now Merge and Process AZFP - glider record



azfp_time = azfp(1).Date;
[~,idx] = unique(azfp_time);
azfp_time(~idx)=azfp_time(~idx)+1e-6;

tgap = 300/86400;
azfp_depth =interp1gap(dcast.time,dcast.depth,azfp_time,tgap);
azfp_pidx = interp1(dcast.time,dcast.profile_index,azfp_time,'nearest');
azfp_lat  = interp1gap(dcast.time,dcast.lat,azfp_time,tgap);
azfp_lon  = interp1gap(dcast.time,dcast.lon,azfp_time,tgap);

ind=find(isnan(azfp_depth));
azfp_depth(ind)=[];
azfp_pidx(ind)=[];
azfp_lat(ind)=[];
azfp_lon(ind)=[];
azfp_time(ind)=[];

echoSv1=azfp(1).Sv;
echoTS1=azfp(1).TS;
range1=azfp(1).TiltCorrRange;range1=range1(1,:);
echoSv1(ind,:)=[];


% sanity check plot AZFP for Profile 61 
% a good since there is data :(
% [xq1,~] = meshgrid(azfp_time,range1);
% [~,yq1] = meshgrid(azfp_depth,range1);
% dq1 = yq1+azfp_depth';
% 
% figure
% pcolor(xq1(:,idx),-dq1(:,idx),echoSv1(idx,:)'); shading flat;
% datetick('x','keeplimits')
% cmocean('balance'); colorbar
% 
% figure; hold on
% plot(dcast.time,-dcast.depth,'k.')
% plot(azfp_time,-azfp_depth,'rx')
% datetick('x','mmm dd','keeplimits') 

%% Next steps
% !$#@ Only about 2-3 days of data :( 
echoSv = echoSv1;
bin_depth = range1;   

%% find indices location of max counts for bottom detection, start at index 10 so it skips the area near transducer
% to eliminate from plots all stuff below seafloor
[~,loc]=max(echoSv(:,10:end),[],2); % why 10?
ExtraY=9; % extra bins to plot beyond bottom % for echo5.mat -> ExtraY=15, otherwise 50
BottomLoc=loc+ExtraY;
[n,l]=max(BottomLoc); % this is to accomodate all different "depths" of profiles <o> in the original one n was BottomLoc(1)

rot_Sv=nan(size(echoSv,1),n);
rot_depth=nan(size(echoSv,1),n);
rot_time=nan(size(echoSv,1),n);

bin_depth=range1;
tmp_depth=azfp_depth*ones(1,n)+ones(length(azfp_depth),1)*bin_depth(1:n);  % dont understand this line.  BUG HERE!! bin_depth(1:n);
tmp_time=repmat(azfp_time,[1,n]);

% % sanity check
subplot(2,1,1)
pcolor(tmp_depth'),shading flat,set(gca,'Ydir','reverse'),colorbar
subplot(2,1,2)
pcolor(echoSv'),shading flat,set(gca,'Ydir','reverse'),colorbar

% starts "rotating" to have the floor signal at the bottom
for ii=1:size(echoSv,1) % first time, second time, so on...
   
    if(BottomLoc(l)-BottomLoc(ii)+1<1) % need to keep this bit otherwise it goes to crap...don't know exactly why
        break;
    end
   
    rot_Sv(ii,BottomLoc(l)-BottomLoc(ii)+1:BottomLoc(l))=echoSv(ii,1:BottomLoc(ii));
    rot_depth(ii,BottomLoc(l)-BottomLoc(ii)+1:BottomLoc(l))=tmp_depth(ii,1:BottomLoc(ii));
    rot_time(ii,BottomLoc(l)-BottomLoc(ii)+1:BottomLoc(l))=tmp_time(ii,1:BottomLoc(ii));
%     
end

clear ii

% x=isnan(rot_Sv);
% rot_Sv(x)=100;
% rot_depth(x)=nan;
%[x y]=ind2sub(size(Output.TS),max_idx); % turns the linear position into i,j index position
%[x y]=ind2sub(size(tmp_depth),ind_d);
% rot_Sv(isnan(rot_Sv))=0;
% rot_depth(isnan(rot_Sv))=0;
% rot_time(isnan(rot_Sv))=0;

figure,plot(rot_depth,'.'),set(gca,'Ydir','reverse'),grid on
figure,set(gcf,'color','w'),pcolor(rot_Sv'),shading flat,set(gca,'Ydir','reverse'),cmocean('balance'),colorbar
% contour(rot_Sv'),set(gca,'Ydir','reverse')
% figure,surf(rot_time,rot_depth,rot_Sv)

% figure
% % pcolor(tmp_time',tmp_depth',rot_Sv'),shading flat
% pcolor(rot_Sv'),shading flat
% % figure,imagesc(rot_Sv')
% % figure,pcolor(rot_depth'),shading flat
% set(gca,'Ydir','reverse')
% cmocean('balance') %balance, curl is greenish
% colorbar
% datetick('x','mmm dd HH PM')
% hold on 
% yyaxis right
% plot(azfp_time,glider_depth,'k.')
% set(gca,'Ydir','reverse')
% datetick('x','mmm dd HH PM')

% %sanity checks
% figure,plot(azfp_time,glider_depth,'.')
% set(gca,'Ydir','reverse')
% datetick('x','HH PM')
% 
% [X,Y]=meshgrid(azfp_time,glider_depth);
% pcolor(X,Y,rot_Sv');shading flat
% contour(rot_time',rot_depth',rot_Sv')
% 
% tmp=rot_depth(1,:);
% 
% [X,Y]=meshgrid(tmp_time,tmp_depth);

%% convert lat/long to decimal degrees ....  why?????
% Lat=((abs(lat)-floor(abs(lat)/100)*100)/60+floor(abs(lat)/100)).*lat./abs(lat);
% Lon=((abs(lon)-floor(abs(lon)/100)*100)/60+floor(abs(lon)/100)).*lon./abs(lon);

%% find the ends of the separate dives and assign each ping to their specific profile
% Kim's solution
% identifies the beginning of each downcast by time difference. Kim changed from 2*N to 25*N to account for variability in ping rate that occassionally occurs
pp=diff(datenum(azfp_time))*24*60*60; % number of seconds between each ping
i_end=find(pp>6); % identify the ends of profiles

% % % using logical indexing instead to find the end of each profile
% tmp=~isnan(glider_depth); % finds all indexes where there are depth values
% pp=diff(tmp); % looks for differences between 0 (a nan) and a 1 (a depth value exists) (x2-x1, i.e. 0-1)
% i_end=(find(pp==-1)); % +1 shows where a downcast starts (transition between a nan and a value will be +), and -1 where it ends (x2-x1)

i_end=[i_end; length(azfp_time)]; % this adds the last index to add the last profile in this file
% i_end(find(diff(i_end)==1))=[]; % this will tell us where each profile ends
profile_nr=[];% allocates a structure for each profile ~4 profiles per hour, it will give a profile number to each ping, i.e. group each ping that belongs in a specific profile

for i=1:length(i_end)
    profile_nr=[profile_nr; ones(i_end(i)-length(profile_nr),1)*i]; %it creates a sctruture with profile nr for each dive
end

clear i

%% combine data from each dive into a profile with .5 m bins
% from glider you need: echo.i_depth=depth, echo.lat, echo.lon each per ping number
% from ProcessAZFP you need:
% bin_depth=range
% Output.Sv and Output.N  -> Rename them to echo.Sv and echo.peak2peak
% Timestamp (per ping)
% echo.i_depth=glider_depth;
% bin_depth=range;
% echo.timestamp=azfp_time;

% clear glider_depth range azfp_time

nr_dives=max(profile_nr);
% dives_depth=floor((min(glider_depth)+min(bin_depth))*2)/2:.3:ceil((max(glider_depth)+max(bin_depth))*2)/2;
%dives_depth=(min(floor((min(rot_depth,[],2))))):.5:(max(ceil(max(rot_depth,[],2))));
dives_depth=0:0.5:600;

% dives_depth=floor((min(glider_depth)+min(bin_depth))*2)/2:.5:ceil((max(glider_depth)+max(bin_depth))*2)/2;
% this makes up a "fake idealized" equally spaced grid that will be use
% to assign echos within a given range; still to find out what the 2 does

% allocating some space...
% dives_lat2=nan(nr_dives,1); % why not just like this????
dives_lat=nan*ones(nr_dives,1);
dives_lon=nan*ones(nr_dives,1);
dives_time=nan*ones(nr_dives,1);
dives_Sv=nan*ones(nr_dives,length(dives_depth));
% dives_nr_averaged=nan*ones(nr_dives,length(dives_depth));

% gets a mean for the entire dive
for i=1:nr_dives %for each dive
    ind_prof=find(profile_nr==i); % find indices that belong to each profile
    dives_time(i)=nanmean(datenum(azfp_time(ind_prof))); % gets a mean time for each dive; can be done without datenum
    dives_lat(i)=nanmean(azfp_lat(ind_prof)); %gets the mean lat/lon for each dive
    dives_lon(i)=nanmean(azfp_lon(ind_prof));
    tmp_Sv=rot_Sv(ind_prof,:);
    tmp_Sv1=echoSv(ind_prof,:); % temporary; gets each echo value for each ping (and all ranges) belonging to each profile/dive (1,2,...)
    % this bit assigns a "true" depth for the bin, taking into consideration to the glider_depth
%     n=size(echoSv,2); %it gets the size of range dimension of azfp bins; it is the same as [~,n]=size(echoSv1)
%     original: temp_dep=echo.i_depth(ind_prof)*ones(1,length(ping{1,1}(nf_rang+1:param.rangebn-1)))+ones(length(ind_prof),1)*bin_depth;
%     tmp_depth=glider_depth(ind_prof)*ones(1,n)+ones(length(ind_prof),1)*bin_depth; %it makes the "real" echo depth matrix by adding the glider depth to the bin range from azfp to get the real depth of each bin
    tmp_depth=rot_depth(ind_prof,:);
    for j=1:length(dives_depth) % the "idealized" depth matrix that will include all ping values within a given range
%         ind_d=find(rot_depth(ind_prof,:)>dives_depth(j)-.25&rot_depth(:)<=dives_depth(j)+.25);
        ind_d=find(tmp_depth(:,:)>dives_depth(j)-.25&tmp_depth(:,:)<=dives_depth(j)+.25); %looks up on the real depth matrix for values within 0.5m from the specific idelized bin
%         dives_Sv(i,j)=10*log10(nanmean(10.^(tmp_Sv(ind_d)/10))); % WHY?????? % takes  the MEAN of ALL VALUES IN THAT DEPTH RANGE IN LINEAR SPACE
        dives_Sv(i,j)=nanmean(tmp_Sv(ind_d));
%         nr_averaged_cells(i,j)=length(find(~isnan(tmp_Sv(ind_d)))); % find number/quantity of averaged cells 
    end
end

figure
set(gcf,'color','w');
h1=imagescn(dives_time,-dives_depth,dives_Sv')
% set(h1,'alphadata',~isnan(dives_Sv'))
datetick('x','mmm dd HH PM','keepticks');
xtickangle(45)%,'Rotation',45.0) 
ylabel('Water depth (m)','FontSize',11)
xlabel('Time UTC (NST-3.5H UTC)','FontSize',11)
axis tight
cmocean('balance')
h2=colorbar; %('AxisLocation','in');
ylabel(h2,'Volume backscatter Sv (dB re m^{-1})','FontSize',11);%,'Rotation',270.0)
set(gca,'YLim',[-600 0])
caxis([-90 -50])