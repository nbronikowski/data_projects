function Z_grid = knnRbfInterpolation(X, Y, Z, gridX, gridY, epsilon, k,regularization)
    % KNNRBFINTERPOLATION Performs k-nearest neighbors interpolation with Radial Basis Functions.
    %
    % X, Y: Coordinates of data points.
    % Z: Values at the data points.
    % gridX, gridY: Coordinates of grid points for interpolation.
    % epsilon: Shape parameter for the RBF, controls the spread.
    % k: The number of nearest neighbors to consider for interpolation.
    % Z_grid: Interpolated values on the grid.
    

    C = regularization;

    % Initialize the output grid
    [a, b]=size(gridX);
    Z_grid = nan(size(gridX));
    
    % Flatten the grid matrices
    gridX = gridX(:);
    gridY = gridY(:);
    
    % Perform the interpolation at each grid point
    for i = 1:length(gridX)
        % Compute distances from the current grid point to all data points
        distances = sqrt((X - gridX(i)).^2 + (Y - gridY(i)).^2);
        
        % Find indices of k nearest neighbors
        [~, indices] = mink(distances, k);
        
        % Retrieve the nearest neighbor coordinates and values
        Xk = X(indices);
        Yk = Y(indices);
        Zk = Z(indices);
        
        % Compute the RBF matrix for the neighbors
        r = sqrt(bsxfun(@minus, Xk, Xk').^2 + bsxfun(@minus, Yk, Yk').^2);
        A = exp(-(epsilon * r).^2);

        % Add a regularization term to the diagonal of the matrix
        A = A + C * eye(size(A));
        
        % Solve the system for weights
        weights = A \ Zk;
        
        % Compute distances from grid point to neighbors
        r_eval = sqrt((gridX(i) - Xk').^2 + (gridY(i) - Yk').^2);
        
        % Evaluate the RBF with the found weights
        Z_grid(i) = exp(-(epsilon * r_eval).^2) * weights;
    end
    
    % Reshape the result back into the original grid shape
    Z_grid = reshape(Z_grid, a,b);
end


