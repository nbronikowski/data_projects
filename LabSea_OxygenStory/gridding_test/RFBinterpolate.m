


function Z_grid=RFBinterpolate(X,Y,Z,gridX,gridY,epsilon)
    % Initialize the output grid
    numGridPoints = numel(gridX);
    Z_grid = zeros(numGridPoints, 1);
    
    % Number of data points
    numDataPoints = length(X);
    
    % Process in chunks to save memory
    chunkSize = 1000;  % Adjust the chunk size based on your system's memory capacity
    numChunks = ceil(numDataPoints / chunkSize);
    
    for i = 1:numChunks
        % Indices for the current chunk
        chunkIndices = (1 + (i-1) * chunkSize):min(i * chunkSize, numDataPoints);
        
        % Distances from grid points to the current chunk of data points
        [chunkGridX, chunkGridY] = ndgrid(gridX, Y(chunkIndices));
        [chunkDataX, chunkDataY] = ndgrid(X(chunkIndices), gridY);
        r = sqrt((chunkGridX - chunkDataX).^2 + (chunkGridY - chunkDataY).^2);
        
        % RBF values for the current chunk
        chunkRBF = exp(-(epsilon * r).^2);
        
        % Partial interpolation result for the current chunk
        Z_grid = Z_grid + sum(chunkRBF .* Z(chunkIndices), 2);
    end
    
    % Reshape to grid
    Z_grid = reshape(Z_grid, size(gridX));
end