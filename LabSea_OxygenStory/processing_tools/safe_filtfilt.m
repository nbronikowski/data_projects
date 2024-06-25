function filtered_data = safe_filtfilt(b, a, data)
    % Find the NaNs
    nanIdx = isnan(data);
    
    % Remove the NaNs
    nonNaNData = data(~nanIdx);
    
    % Apply filtfilt
    nonNaNFiltered = filtfilt(b, a, nonNaNData);
    
%     nonNaNFiltered =medfilt1(nonNaNData,3);
    
    % Restore the NaNs
    filtered_data = nan(size(data));
    filtered_data(~nanIdx) = nonNaNFiltered;
end
