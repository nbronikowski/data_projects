function [subsetData,cols_to_keep]=deleteAlmostEmptyColumns(dat,depthg)

[~, N] = size(dat);
% Initialize a logical vector to identify columns to keep
cols_to_keep = false(1, N);

% Loop through each column to find the data span
for t = 1:N
    non_nan_idx = find(~isnan(dat(:, t))); % Find non-NaN indices
    
    % Check if there are any non-NaN entries
    if ~isempty(non_nan_idx)
        data_span = max(depthg(non_nan_idx)) - min(depthg(non_nan_idx)) + 1; % Compute data span
        
        % Check if data span is at least 50
        if data_span >= 50
            cols_to_keep(t) = true;
        end
    end
end

% Keep only columns with sufficient data span
subsetData = dat(:, cols_to_keep);
idx = cols_to_keep;