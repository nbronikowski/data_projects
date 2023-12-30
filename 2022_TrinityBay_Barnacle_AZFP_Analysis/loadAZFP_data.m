clear; close all; clc

azfp_tbox_path = ['/Users/research/Library/' ...
    'Mobile Documents/com~apple~CloudDocs/toolbox/AZFP_MATLAB_Toolbox_v1.9/'];
addpath(azfp_tbox_path);
Parameters.ProcDir = 0;
% Parameters.NoiseFloor = [0 0 0 0];
Parameters.Salinity = 30; 
Parameters.Bins2Avg = 15; % approx 0.5 m range
Parameters.UseTiltCorr = 1;
Parameters.Time2Avg = 1;
Parameters.Pressure = 50;
Parameters.Plot = 0;
Parameters.Value2Plot = 2;


azfp_flist_dir = './raw_azfp_files_20221120/';
flist = dir([azfp_flist_dir, '*.01*']);

xml_flist_dir = './raw_azfp_files_20221120/xmlfiles/';
xml_list = dir([xml_flist_dir, '*.XML']);

% Get the names of XML files
xml_fnames = {xml_list.name};
xml_fnames = sort(xml_fnames);

var_list = [
    {'Date'                }
    {'Tx'                  }
    {'Ty'                  }
    {'T'                   }
    {'BatteryMain'         }
    {'BatteryTx'           }
    {'Depth'               }
    {'filename'            }
    {'HourlyAvgTemp'       }
    {'SoundSpeed'          }
    {'N'                   }
    {'Range'               }
    {'TiltCorrRange'       }
    {'Sv'                  }
    {'TS'                  }
    {'seaAbs'              }
    {'Freq'                }
    {'dirname'             }
    {'StartEndRangeNumBins'}
    {'Bins2Avg'            }
    {'Time2Avg'            }
    {'BurstInt'            }
    {'PingPerProfile'      }
    {'NumAcqPings'         }
    {'DataType'            }
    {'PulseLength'         }
    {'BinInc'              }
    {'TimeInc'             }
    {'SerialNumber'        }];
for jj = 1:length(var_list)
    for ll=1:4
        azfp(ll).(var_list{jj})=[];
    end
end

% Load Glider CTD Data
load('glider_ctd_data.mat');
% dat.profile_index(mod(dat.profile_index(~isnan(dat.profile_index)),1)==0)=NaN;


% Set up a default XML file (first in the list)
defaultXML = xml_fnames{5};
Parameters.xmlpathname = fullfile(xml_list(1).folder);
Parameters.datapathname = fullfile(flist(1).folder);
for i = 1:length(flist)
    % Extract the date part of the filename
    fileDate = flist(i).name(1:6);

    % Initialize the matchedXML variable
    matchedXML = defaultXML;

    % Iterate over XML filenames to find a match
    for j = 1:length(xml_fnames)
        xmlDate = xml_fnames{j}(1:6);
        if strcmp(fileDate, xmlDate)
            matchedXML = xml_fnames{j};
            break; % Exit the loop once a match is found
        end
    end

    % Use the matched or default XML filename
    Parameters.datafilename =fullfile(flist(i).name);
    Parameters.xmlfilename = fullfile(matchedXML);
    [Output, Par, PavgArr] = ProcessAZFP(Parameters);

%     t1 = min(Output(1).Date); t2 = max(Output(1).Date);
%     
%     disp([datestr(t1),' - ',datestr(t2)])
%     gl_idx = dat.time>t1 & dat.time<t2;
% 
%     figure(1); hold on
%     plot(dat.time,-dat.depth)
%     plot(dat.time(gl_idx),-dat.depth(gl_idx),'.r','LineWidth',2)
%     xlim([t1-0.1 t2+0.1])
%     close all
    
    for jj = 1:length(var_list)
        for ll = 1:length(Output)
            temp = Output(ll).(var_list{jj}); 
            if i == 1
                azfp(ll).(var_list{jj})=[azfp(ll).(var_list{jj}),temp];
            else
                azfp(ll).(var_list{jj})=cat(1,azfp(ll).(var_list{jj}),temp);
            end
        end
    end
end

save('azfp_out.mat','azfp','-v7.3')