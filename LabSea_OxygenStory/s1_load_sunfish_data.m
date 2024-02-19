clear; close all; clc;

% Set ERDDAP base URL
baseURL = 'https://www.smartatlantic.ca/erddap/tabledap/';

datasetID = 'sunfish_labrador_sea_2022_60ac_212a_97de';
dataFormat = '.mat'; % The format we want to download

% Construct the query URL for your specific data request
% List of variables to include in the query
variables = [...
    'time,time_ctd,latitude,longitude,wpt_lat,wpt_lon,' ...
    'pitch,roll,heading,water_u,water_v,profile_index,' ...
    'profile_direction,depth,temperature,pressure,conductivity,salinity,' ...
    'abs_salinity,density,oxygen_calphase,' ...
    'cdom_ref,chlor_ref,bb700_ref,cdom_sig,chlor_sig,bb700_sig'...
            ];

% Time constraints for the query
startTime = '2021-12-05T00:00:00Z'; % start of 53 N LabSea data
endTime = '2022-06-8T00:00:00Z';   % end of 53 N LabSea data

options = weboptions('Timeout', 300); % Increase timeout to 5 min

% Encode the variables for the URL (replace spaces with %2C for commas)
encodedVars = strrep(variables, ' ', '%2C');

% Construct the full URL for the data request
fullURL = sprintf('%s%s%s?%s&time>=%s&time<=%s', ...
                  baseURL, datasetID, dataFormat, encodedVars, startTime, endTime);

% Specify the local filename to save the downloaded file
localFilename = 'sunfish_data.mat'; % The filename to save the downloaded .mat file
% Check if the 'mat_files' folder exists, if not, create it
folderName = 'mat_files'; % Define the folder name
if ~exist(folderName, 'dir')
    mkdir(folderName); % Create the folder if it doesn't exist
end

% Specify the local filename to save the downloaded file, including folder path
localFilename = fullfile(folderName,localFilename); % Save inside the 'mat_files' folder

% Use websave to download the file
try
    websave(localFilename, fullURL,options);
    fprintf('Data successfully downloaded to %s\n', localFilename);
catch webException
    fprintf('Error downloading data: %s\n', webException.message);
    return;
end

% Now that you've downloaded the file, read it into MATLAB
% Do a bit of cleaning to make both data sets similar
load(localFilename)
sunfish = sunfish_labrador_sea_2022_60ac_;
clear sunfish_labrador_sea_2022_60ac_
sunfish.conductivity = sunfish.conductivity./10;
sunfish.u = sunfish.water_u;
sunfish.v = sunfish.water_v;
sunfish=rmfield(sunfish,'water_u');
sunfish=rmfield(sunfish,'water_v');
save(localFilename,'sunfish')




