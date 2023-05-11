% Clear the workspace and command window
clc; clear; close all;

try
    % Open the CSV file and read in all the lines
    fid = fopen("Results March24 2023.csv");
    Llines = {};
    tline = fgetl(fid);
    while ischar(tline)
        Llines{end+1,1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);

    % Define the delimiter used in the CSV file
    delimiter = ',';

    % Split each row by comma delimiter 
    func = @(input)strsplit(input,delimiter);
    % Apply this function to each row in the cell array to get the final data
    % As each row is of different length, we set the UniformOutput to false
    raw_data = cellfun(func,Llines,'UniformOutput',false);

    % Get the column headers
    var_names = raw_data{1};
    N = length(var_names);
    M = length(Llines);

    % Find the row index for the column headers
    header_idx = NaN(M,1);
    for i = 1:M
        cl_header = string(raw_data{i,:});
        var_query = string({'File #'});
        if ~isempty(find(strcmp(cl_header,var_query)))
            header_idx(i)=i;
        end
    end
    header_idx(isnan(header_idx))=[];

    % Preallocate cell array for sorted data
    sorteD = cell(M-length(header_idx),N);

    % Loop through each header and its associated data, and sort them into the correct position
    t = 1;
    for j = 1:length(header_idx)
        header_row = raw_data{header_idx(j)};
        data_start_row = header_idx(j)+1; % Index of the first row of data after the header
        data_end_row = M; % Index of the last row of data
        
        % Check if there is another header row in the data
        if j < length(header_idx)
            data_end_row = header_idx(j+1)-1; % Index of the last row of data before the next header
        end
        
        % Loop through each row of data and sort it based on the header
        for k = data_start_row:data_end_row
            current_row_data = raw_data{k};
            for s = 1:length(header_row)
                var_query = string(header_row{s});
                var_list = string(var_names);
                master_header_position = find(strcmp(var_list,var_query));
                if ~isempty(master_header_position)
                    sorteD(t, master_header_position) = current_row_data(s);
                end
            end
            t = t+1;
        end
    end

    % Write the sorted data to a new CSV file
    writetable(cell2table(sorteD,'