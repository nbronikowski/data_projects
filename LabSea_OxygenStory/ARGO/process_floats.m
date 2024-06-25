clear; clc; close all;

fpath = './new_nc/';
flist=dir([fpath,'*.nc']);

out = ncinfo(fullfile([fpath,flist(40).name]));

VAR_LIST = ...
{'TIME',... % 'days since 1950-01-01T00:00:00Z'
'LATITUDE',...
'LONGITUDE',...
'PRES',...
'PRES_QC',...
'PRES_ADJUSTED',...
'PRES_ADJUSTED_QC',...
'PSAL',...
'PSAL_QC',...
'PSAL_ADJUSTED',...
'PSAL_ADJUSTED_QC',...
'TEMP',...
'TEMP_QC',...
'TEMP_ADJUSTED',...
'TEMP_ADJUSTED_QC',...
'DOX2_ADJUSTED',...
'DOX2_ADJUSTED_QC'};

% check if float has oxygen data and then copy data to a MAT Structure
PG = 0:1:1029; PG=PG(:);

float_map = struct();
for i = 1:length(VAR_LIST)
    float_map.(VAR_LIST{i})=[];
end
float_map.WMO_ID=[];

sfloat = [];
tfloat = [];
for i = 1:length(flist)
    out = ncinfo(fullfile([fpath,flist(i).name]));
    if any(strcmp({out.Variables.Name},'DOX2_ADJUSTED'))
        %% GRAB STUFF FROM THIS FLOAT NEAR THE GLIDER
        for j = 1:length(VAR_LIST)
            if any(strcmp({out.Variables.Name},VAR_LIST{j}))
                tfloat.(VAR_LIST{j}) = double(ncread(fullfile([fpath,flist(i).name]),VAR_LIST{j}));
            else
                tfloat.(VAR_LIST{j}) = NaN*double(ncread(fullfile([fpath,flist(i).name]),'DOX2_ADJUSTED'));
            end
        end

        %% Interpolate to common grid
        for j = 1:length(VAR_LIST)
            if strcmp(VAR_LIST{j},'TIME') || strcmp(VAR_LIST{j},'LONGITUDE') || strcmp(VAR_LIST{j},'LATITUDE')
                % Get the total number of elements in the array
                M = numel(tfloat.(VAR_LIST{j}));
                % Reshape the array into 1xM format
                sfloat.(VAR_LIST{j}) = reshape(tfloat.(VAR_LIST{j}), 1,M);
            else
                flags = find(~isnan(squeeze(tfloat.(VAR_LIST{j}(:)))));
                if length(flags)>20
                    sfloat.(VAR_LIST{j}) = gridcolumns(tfloat.(VAR_LIST{j}),tfloat.PRES,tfloat.PRES_QC,PG,1);
                else
                    sfloat.(VAR_LIST{j}) = NaN(length(PG),length(tfloat.TIME));
                end
            end
        end
        sfloat.TIME = sfloat.TIME + datenum(1950,1,1);
        sfloat.WMO_ID = zeros(size(sfloat.TEMP))+str2double(ncreadatt(fullfile([fpath,flist(i).name]),'/','wmo_platform_code'));
        
        % ADD RESULTS
        for j = 1:length(VAR_LIST)
            float_map.(VAR_LIST{j}) = [float_map.(VAR_LIST{j}),sfloat.(VAR_LIST{j})];
        end
        float_map.WMO_ID = [float_map.WMO_ID,sfloat.WMO_ID];
        sfloat = []; tfloat = [];
    end    
end

%% delete empty columns

[~,colIdx]=deleteAlmostEmptyColumns(float_map.DOX2_ADJUSTED,PG);
VAR_LIST = fieldnames(float_map);
float_map2=struct();
for j = 1:length(VAR_LIST)
    float_map2.(VAR_LIST{j})=float_map.(VAR_LIST{j})(:,colIdx);
end
float_map = float_map2; clear float_map2

figure()
ti=tiledlayout(3,1);

nexttile
T = float_map.TEMP_ADJUSTED;
T(float_map.TEMP_ADJUSTED_QC>1)=NaN;
imagescn(float_map.TIME,-PG,T)

nexttile
S = float_map.PSAL_ADJUSTED;
S(float_map.PSAL_ADJUSTED_QC>1)=NaN;
imagescn(float_map.TIME,-PG,S);

nexttile
O2 = float_map.DOX2_ADJUSTED;
O2(float_map.DOX2_ADJUSTED_QC>1)=NaN;
imagescn(float_map.TIME,-PG,O2);


float_map.DOX2_ADJUSTED(float_map.DOX2_ADJUSTED_QC>1)=NaN;

%% plot
ufloats = unique(float_map.WMO_ID(:));
colList = seminfhaxby(length(ufloats)+10);

clear h
figure()
t=tiledlayout(1,1);
nexttile; hold on
topoplot([-65 -44 45 70],[],[0:-1/2:-5],'-k')
for i = 1:length(ufloats)
   id = float_map.WMO_ID(1,:)==ufloats(i);
   h(i)=plot(float_map.LONGITUDE(id),float_map.LATITUDE(id),'.','MarkerSize',10,'Color',colList(i,:)); 
   legendInfo{i}=num2str(ufloats(i));
end
legend(h,legendInfo,'location','bestoutside')
formatplot
t.TileSpacing='compact';
subtitle(t,'Map of Argo Oxygen Profiles 2020 - 2022')
save_figure(gcf,'ArgoFloatMap',[6 7],'.png','300')


%% Save MAT file
argo = float_map;
argo.pgrid = PG;
save('argo_oxygen.mat','argo');



