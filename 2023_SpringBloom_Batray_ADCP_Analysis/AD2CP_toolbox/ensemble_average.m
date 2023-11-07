function ds_avg = ensemble_average(ds, secs)
    field_names = fieldnames(ds);
    ds_avg = struct();
    num_chunks = floor((ds.time(end) - ds.time(1))*86400 / secs);
    
    % Initialize a vector to hold the averaged time data
    avg_time = zeros(num_chunks, 1);
    
    % Loop through each chunk and calculate the mean time
    for j = 1:num_chunks
        idx_start = (j-1)*secs + 1;
        idx_end = j*secs;
        if idx_end > length(ds.time)
            idx_end = length(ds.time);
        end
        chunk_time = ds.time(idx_start:idx_end);
        avg_time(j) = nanmean(chunk_time);
    end
    
    % Assign the averaged time data to the new structure
    ds_avg.time = avg_time;
    
    % Loop through each field in the structure
    for i = 1:length(field_names)
        field_name = field_names{i};
        
        % Skip the time field since we've already handled it
        if strcmp(field_name, 'time')
            continue;
        end
        
        data = ds.(field_name);
        
        % Ensure data has a time dimension
        assert(size(data,2) >= 1, 'Data must have a time dimension');
        
        % Initialize a matrix to hold the averaged data
        avg_data = zeros(num_chunks, size(data, 1));
        
        % Loop through each chunk and calculate the mean
        for j = 1:num_chunks
            idx_start = (j-1)*secs + 1;
            idx_end = j*secs;
            if idx_end > size(data, 2)
                idx_end = size(data, 2);
            end
            chunk_data = data(:, idx_start:idx_end);
            avg_data(j, :) = nanmean(chunk_data, 2);
        end
        
        % Assign the averaged data to the new structure
        ds_avg.(field_name) = avg_data';
    end
end