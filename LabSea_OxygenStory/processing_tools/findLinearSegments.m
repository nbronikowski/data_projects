function linearIndices = findLinearSegments(data, tolerance)

    % Initialize variables
    n = length(data);
    linearIndices = []; % Vector to store indices of linear segments
    
    % Calculate slopes (differences) between consecutive points
    slopes = [diff(data); 0]; % Append 0 to make slopes array same length as data
    
    % Iterate through data to find indices part of linear segments
    for i = 2:n-1 % Skip the first and last data points for comparison
        if isnan(data(i-1)) || isnan(data(i)) || isnan(data(i+1))
            continue; % Skip if current or neighboring data points are NaN
        end
        
        % Check for consistent slope within tolerance
        if abs(slopes(i-1) - slopes(i)) <= tolerance && abs(slopes(i) - slopes(i+1)) <= tolerance
            linearIndices = [linearIndices; i]; % Add index if part of linear segment
        end
    end
    
    % Include start and end points of linear segments by checking neighbors
    if ~isempty(linearIndices)
        % Check if the second point in a segment is linear and include the first one
        if linearIndices(1) > 1 && ~ismember(linearIndices(1)-1, linearIndices)
            if abs(slopes(linearIndices(1)-1) - slopes(linearIndices(1))) <= tolerance
                linearIndices = [linearIndices(1)-1; linearIndices];
            end
        end
        
        % Check if the second-to-last point in a segment is linear and include the last one
        if linearIndices(end) < n && ~ismember(linearIndices(end)+1, linearIndices)
            if abs(slopes(linearIndices(end)) - slopes(linearIndices(end)+1)) <= tolerance
                linearIndices = [linearIndices; linearIndices(end)+1];
            end
        end
    end
% % Initialize variables
%     n = length(data);
%     isLinear = false(n-1, 1); % Logical array to mark linear differences
%     
%     % Calculate slopes (differences) between consecutive points, ignoring NaNs
%     slopes = diff(data);
%     
%     % Iterate through slopes to find linear segments
%     for i = 1:n-2
%         if isnan(slopes(i)) || isnan(slopes(i+1))
%             continue; % Skip NaNs
%         end
%         if abs(slopes(i) - slopes(i+1)) <= tolerance
%             isLinear(i) = true; % Mark as part of a linear segment
%         end
%     end
%     
%     % Identify starting and ending indices of linear segments
%     linearSegments = [];
%     i = 1;
%     while i <= length(isLinear)
%         if isLinear(i)
%             startIdx = i;
%             while i <= length(isLinear) && isLinear(i)
%                 i = i + 1;
%             end
%             endIdx = i; % The end of the linear segment
%             linearSegments = [linearSegments; [startIdx, endIdx]]; % Append the indices
%         else
%             i = i + 1;
%         end
%     end
%     
%     % Adjust indices to match original data array
%     linearSegments = linearSegments + 1;
%     if ~isempty(linearSegments)
%         linearSegments(:, 2) = linearSegments(:, 2) + 1; % Adjust end index for inclusive range
%     end
end
