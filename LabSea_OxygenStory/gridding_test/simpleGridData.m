
function Z_grid = simpleGridData(X, Y, Z, gridX,gridY)
    % Initialize the output grid with NaNs
%     Z_grid = nan(size(gridX));
%     
%     % Calculate gradients in each dimension
%     F = scatteredInterpolant(X, Y, Z, 'linear', 'none');
%     Z_grid = F(gridX, gridY);
    Z_grid = griddata(X,Y,Z,gridX,gridY,'nearest');

end