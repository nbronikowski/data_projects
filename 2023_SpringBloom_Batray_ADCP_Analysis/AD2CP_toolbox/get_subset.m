function ds_subset=get_subset(ds,t1,t2)  
    fn = fieldnames(ds);
    idx = find(ds.time>=t1 & ds.time<=t2);
    ds_subset = struct();
    for i = 1:length(fn)
        if size(ds.(fn{i}), 2) == length(ds.time)  % Check if time is the second dimension
            ds_subset.(fn{i}) = ds.(fn{i})(:, idx);  % Select only the columns corresponding to the desired times
        elseif isequal(ds.(fn{i}), ds.time)  % Check if the field is the time vector itself
            ds_subset.(fn{i}) = ds.(fn{i})(idx);  % Select only the desired times
        else
            warning('Field %s does not match the time dimension and will not be subsetted.', fn{i});
        end
    end
end