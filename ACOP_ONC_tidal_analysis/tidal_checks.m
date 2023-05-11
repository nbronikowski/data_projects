clear; clc; close all;

flist = dir('./search28221436/*.mat');
load([flist.folder,'/',flist.name]);

elev=data(11).dat;
elev=elev-nanmean(elev);
time_x=data(11).time;
tintval = nanmedian(diff(time_x))*24; % 120 seconds


addpath(genpath("t_tide_v1.4beta/"));

% infername=['P1';'K2'];
% inferfrom=['K1';'S2'];
% infamp=[.33093;.27215];
% infphase=[-7.07;-22.40];

[tidestruc,pout]=t_tide(elev,...
    'interval',tintval, ...                     % hourly data
    'start',time_x(1),...                       % start time is datestr(tuk_time(1))
    'latitude',metadata.location.lat_degrees);  % Latitude of obs

time_x = (time_x-time_x(1)); % time in days


clf;orient tall;
set(gcf,'defaultaxestickdir','out','defaultaxestickdirmode','manual');
subplot(411);
plot(time_x,elev,'b');
line(time_x,pout,'color',[0 .5 0]);
line(time_x,elev-pout,'linewi',2,'color','r');
xlabel('Days in 2022');
ylabel('Elevation (m)');
legend('Original Time series','Tidal prediction from Analysis',...
    'Original time series minus Prediction','location','best');
title('Tidal analysis with T\_tide (Holyrood Subsea Observatory)');

ax(1)=subplot(412);
fsig=tidestruc.tidecon(:,1)>tidestruc.tidecon(:,2); % Significant peaks
semilogy([tidestruc.freq(~fsig),tidestruc.freq(~fsig)]',[.0005*ones(sum(~fsig),1),tidestruc.tidecon(~fsig,1)]','.-r');
line([tidestruc.freq(fsig),tidestruc.freq(fsig)]',[.0005*ones(sum(fsig),1),tidestruc.tidecon(fsig,1)]','marker','.','color','b');
line(tidestruc.freq,tidestruc.tidecon(:,2),'linestyle',':','color',[0 .5 0]);
set(gca,'ylim',[.0005 1],'xlim',[0 .5]);
xlabel('frequency (cycles/hour)');
text(tidestruc.freq,tidestruc.tidecon(:,1),tidestruc.name,'rotation',45,'vertical','base');
ylabel('Amplitude (m)');
text(.3,.4,'Analyzed lines with 95% significance level','fontweight','bold');
text(.3,.2,'Significant Constituents','color','b');
text(.3,.1,'Insignificant Constituents','color','r');
text(.3,.05,'95% Significance Level','color',[0 .5 0]);

ax(2)=subplot(413);
errorbar(tidestruc.freq(~fsig),tidestruc.tidecon(~fsig,3),tidestruc.tidecon(~fsig,4),'.r');
hold on;
errorbar(tidestruc.freq(fsig),tidestruc.tidecon(fsig,3),tidestruc.tidecon(fsig,4),'o');
hold off;
set(gca,'ylim',[-45 360+45],'xlim',[0 .5],'ytick',[0:90:360]);
xlabel('frequency (cycles/hour)');
ylabel('Greenwich Phase (deg)');
text(.3,330,'Analyzed Phase angles with 95% CI','fontweight','bold');
text(.3,290,'Significant Constituents','color','b');
text(.3,250,'Insignificant Constituents','color','r');

ax(3)=subplot(414);
ysig=elev;
yerr=elev-pout;
nfft=389;
bd=isnan(ysig);
gd=find(~bd);
bd([1:(min(gd)-1) (max(gd)+1):end])=0;
ysig(bd)=interp1(gd,ysig(gd),find(bd)); 
%[Pxs,F]=psd(ysig(isfinite(ysig)),nfft,1,[],ceil(nfft/2));
[Pxs,F]=pwelch(ysig(isfinite(ysig)),hanning(nfft),ceil(nfft/2),nfft,1);
Pxs=Pxs/2;
%%[Pxso,Fo]=psd(ysig(isfinite(ysig)),nfft,1,[],ceil(nfft/2));

%[Pxs,F]=pmtm(ysig(isfinite(ysig)),4,4096,1);
yerr(bd)=interp1(gd,yerr(gd),find(bd)); 
%[Pxe,F]=psd(yerr(isfinite(ysig)),nfft,1,[],ceil(nfft/2));
[Pxe,F]=pwelch(yerr(isfinite(ysig)),hanning(nfft),ceil(nfft/2),nfft,1);
Pxe=Pxe/2;
%[Pxe,F]=pmtm(yerr(isfinite(ysig)),4,4096,1);

semilogy(F,Pxs);
line(F,Pxe,'color','r');
xlabel('frequency (cph)');
ylabel('m^2/cph');
% text(.2,1e2,'Spectral Estimates before and after removal of tidal energy','fontweight','bold');
% text(.3,1e1,'Original (interpolated) series','color','b');
% text(.3,1e0,'Analyzed Non-tidal Energy','color','r');

% Zoom  on phase/amplitude together
linkaxes(ax,'x');



