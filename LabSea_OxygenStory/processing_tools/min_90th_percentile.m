function min_90th = min_90th_percentile(data)
    % Sort the data in ascending order, excluding NaN values
    sorted_data = sort(data(~isnan(data)));

    % Calculate the index corresponding to the 90th percentile
    percentile_index = ceil(0.90 * length(sorted_data));

    % Calculate the 90th percentile value
    percentile_90 = sorted_data(percentile_index);

    % Filter the data to include only values less than or equal to the 90th percentile
    filtered_data_90 = data(data <= percentile_90);

    % Calculate the minimum value of the 90% of the data
    min_90th = nanmin(filtered_data_90);
end