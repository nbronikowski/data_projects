function min_75th = min_75th_percentile(data)
    % Sort the data in ascending order, excluding NaN values
    sorted_data = sort(data(~isnan(data)));

    % Calculate the index corresponding to the 75th percentile
    percentile_index = ceil(0.75 * length(sorted_data));

    % Calculate the 75th percentile value
    percentile_75 = sorted_data(percentile_index);

    % Filter the data to include only values less than or equal to the 75th percentile
    filtered_data_75 = data(data <= percentile_75);

    % Calculate the minimum value of the 75% of the data
    min_75th = nanmin(filtered_data_75);
end