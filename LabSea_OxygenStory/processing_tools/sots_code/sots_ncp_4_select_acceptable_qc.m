% Subsamples the mooring data over the longest possible continuous time 
% period of qc values less than qc_limit

%%

% selects all the variables in the structure 
var_names = fieldnames(mooring_data);

% creates an empty cell array to store the variables that have associated
% qc vectors
var_names_check_qc = {};

% creates an empty vector to store the qc vectors of each of the variables
var_names_qc_values = [];

% populates var_names_check_qc
for i = 1:length(fieldnames(mooring_data))
    
    if any(contains(var_names,strcat(var_names{i},'_qc'))) 
        
        var_names_check_qc{end+1} = var_names{i};
        
    end
    
end

% Sets the initial values for the start and end of acceptable data at the start and end of the timeseries
start_good = 1;

end_good = length(mooring_data.time);

% Loops through each of the variables remaining
for i = 1:length(var_names_check_qc)
    
        % calculates a logical vector of qc data below qc_limit, then calculates
        % stepwise changes in this vector. Accordingly, changes of 1 will 
        % indicate data become acceptable, and changes of -1 will indicate
        % data becoming unacceptable
        qc_change_idx = diff(mooring_data.(strcat(var_names_check_qc{i},'_qc'))<qc_limit);
    
        % sets start_good to the highest of: itself or the highest index 
        % of a change of 1 for the current variable
        if~(any(qc_change_idx==1)) 

            start_good = 1;

        else

            start_good = max(start_good,find(qc_change_idx==1,1));

        end

        
        % sets end_good to the lowest of: itself or the lowest index 
        % of a change of 1 for the current variable
        end_good = min(end_good,find(qc_change_idx==-1,1));
        
        % concatenates the vector of qc values of the current variable to var_names_qc_values
        var_names_qc_values = [var_names_qc_values mooring_data.(strcat(var_names_check_qc{i},'_qc'))];

end

% trims var_names_qc_values to (start_good:end_good)
var_names_qc_values = var_names_qc_values(start_good:end_good);

% loops through each variable in mooring_data
for i = 1:length(var_names)
    
    % trims their length in the structure to (start_good:end_good)
    mooring_data.(var_names{i}) = mooring_data.(var_names{i})(start_good:end_good);
    
end

disp(['Acceptable data starts at ',datestr(mooring_data.time(1))])

disp(['Acceptable data end at ',datestr(mooring_data.time(end))])
