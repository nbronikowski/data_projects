
load('glider_data.mat')
% 
% tstamp = gdat.timeDateNum;
% values = [gdat.cond*10,gdat.temp,gdat.pres,...
%           gdat.oxygen_calphase,gdat.oxygen_temp,gdat.oxygen_concentration];
% idnan = find(isnan(tstamp) | isnan(gdat.pres) | isnan(gdat.cond));
% tstamp(idnan)=[];
% values(idnan,:)=[];
% 
% channel = {'Conductivity','Temperature','Sea Pressure','O2 Phase','O2 Temp','O2 Concentration'};
% unit = {'mS/cm','C','dbar','degrees','C','umol/L'};
% rsk = RSKcreate('tstamp',tstamp,'values',values,...
%                'channel',channel,'unit',unit,'model','glider');
% [rsk,hasProfile] = RSKfindprofiles(rsk);
% rsk = RSKtimeseries2profiles(rsk);
% 
% rsk = RSKderivedepth(rsk,'latitude',mean(gdat.lat,'omitnan'));
% rsk = RSKderivevelocity(rsk);
% rsk = RSKremoveloops(rsk,'threshold',0.05);
% rsk = RSKderivesalinity(rsk);
% 
% % lag = RSKcalculateCTlag(rsk);
% rsk = RSKalignchannel(rsk,'channel','temperature','lag',-1);
% rsk = RSKderivesalinity(rsk);
% 
% RSKplotprofiles(rsk,'channel',{'Salinity'},'Profile',[3 4])
% 
% [rsk,~]=RSKdespike(rsk,'channel','salinity','windowlength',7);
% 
% 
% %% Store Results Salinity Correction
% channel = {rsk.channels.longName};
% profidx = rsk.profiles.originalindex;
% channel_o = strrep(channel,' ','');
% temp_out.time = [];
% for j = 1:length(channel_o)
%     temp_out.(channel_o{j})=[];
% end
% 
% for i = 1:length(profidx)
%     temp_out.time = [temp_out.time;rsk.data(profidx(i)).tstamp];
%     for j = 1:length(channel_o)
%         temp_vals = rsk.data(profidx(i)).values(:,j);
%         temp_out.(channel_o{j}) = [temp_out.(channel_o{j});temp_vals];
%     end
% end
% % MERGE with dat for Salinity
% maxgapval = nanmedian(diff(temp_out.time))*20;
% [~,id] = unique(temp_out.time);
% gdat.salt_corrected = interp1gap(temp_out.time(id),temp_out.Salinity(id),gdat.timeDateNum,maxgapval);
% % gdat.temp = interp1gap(temp_out.time(id),temp_out.Temperature(id),gdat.timeDateNum,maxgapval);
% gdat.y_velocity = interp1gap(temp_out.time(id),temp_out.Velocity(id),gdat.timeDateNum,maxgapval);
% 
% save('./../glider_data.mat','gdat','-v7.3');

idnan = ~(isnan(gdat.temp) | isnan(gdat.cond));
tempNaN = gdat.temp(idnan);
saltNaN = gdat.salt(idnan);
condNaN = gdat.cond(idnan);
presNaN = gdat.pres(idnan);
minAdv = -2.5; maxAdv = 2.5;
[bestAdv,bestSal,advances,goodness,binIndex,salinities] = find_cond_advance(tempNaN,condNaN,presNaN,...
    minAdv,maxAdv,'showprogress','yes');

gdat.salt_corrected = NaN*gdat.salt;
gdat.salt_corrected(idnan)=bestSal;

save('./../glider_data.mat','gdat','-v7.3');

