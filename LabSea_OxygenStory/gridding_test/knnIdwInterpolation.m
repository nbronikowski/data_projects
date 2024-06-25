function Z_grid = knnIdwInterpolation(X, Y, Z, gridXo, gridYo, k)
    % KNNIDWINTERPOLATION Performs k-nearest neighbors interpolation with inverse distance weighting.
    %
    % X, Y: Coordinates of data points.
    % Z: Values at the data points.
    % gridX, gridY: Coordinates of grid points for interpolation.
    % k: The number of nearest neighbors to consider for interpolation.
    % Z_grid: Interpolated values on the grid.

    % Initialize the output grid
    Z_grid = nan(size(gridXo));
    
    % Flatten the grid matrices
    gridX = gridXo(:);
    gridY = gridYo(:);
    
    % Preallocate space for distances and indices
    distances = zeros(length(X), 1);
    indices = zeros(length(X), 1);
    
    % Perform the interpolation at each grid point
    for i = 1:length(gridX)
        % Compute distances from the current grid point to all data points
        for j = 1:length(X)
            distances(j) = sqrt((X(j) - gridX(i))^2 + (Y(j) - gridY(i))^2);
        end
        
        % Find indices of k nearest neighbors
        [~, indices] = mink(distances, k);
        
        % Retrieve the nearest neighbor distances and values
        nearestDists = distances(indices);
        nearestValues = Z(indices);
        
        % Compute weights using inverse distance weighting
        weights = 1 ./ (nearestDists + eps); % Add a small value to avoid division by zero
        
        % Compute weighted average
        Z_grid(i) = sum(weights .* nearestValues) / sum(weights);
    end
    
    % Reshape the result back into the original grid shape
    Z_grid = reshape(Z_grid, size(gridXo));
end