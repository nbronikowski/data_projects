clear; close all; clc;

% Set ERDDAP base URL
baseURL = 'https://www.smartatlantic.ca/erddap/tabledap/';
datasetID = 'mun_glider_data_pearldiver_labrador_sea_2019';
dataFormat = '.mat'; % The format we want to download

% Construct the query URL for your specific data request
% List of variables to include in the query
variables = [...
    'time,time_ctd,latitude,longitude,wpt_lat,wpt_lon,' ...
    'pitch,roll,heading,u,v,profile_index,' ...
    'profile_direction,depth,temperature,pressure,conductivity,salinity,' ...
    'abs_salinity,density,oxygen_calphase'];

% Time constraints for the query
startTime = '2020-01-15T00:00:00Z'; % start of Central LabSea data
endTime = '2020-05-25T00:00:00Z';   % end of Central LabSea data

options = weboptions('Timeout', 300); % Increase timeout to 5 min

% Encode the variables for the URL (replace spaces with %2C for commas)
encodedVars = strrep(variables, ' ', '%2C');

% Construct the full URL for the data request
fullURL = sprintf('%s%s%s?%s&time>=%s&time<=%s', ...
                  baseURL, datasetID, dataFormat, encodedVars, startTime, endTime);

% Specify the local filename to save the downloaded file
localFilename = 'pearldiver_data.mat'; % The filename to save the downloaded .mat file
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
load(localFilename)

pearldiver = mun_glider_data_pearldiver_labr;
clear mun_glider_data_pearldiver_labr
save(localFilename,'pearldiver')