function [vel_referenced, bin_centers, vel_referenced_std] = shear_method(U, V, W, vx, vy, bins, depth, dz)
    % Filter out NaN rows and columns
    nanind_rows = all(isnan(U), 2);
    U(nanind_rows, :) = [];
    V(nanind_rows, :) = [];
    W(nanind_rows, :) = [];
    bins(nanind_rows) = [];
    
    nanind_cols = all(isnan(U), 1);
    U(:, nanind_cols) = [];
    V(:, nanind_cols) = [];
    W(:, nanind_cols) = [];
    depth(nanind_cols) = [];
    
    % Calculate bin size and number
    bin_size = diff(bins(1:2));
    bin_num = length(bins);
    
    % Calculate actual depths of the ADCP bins
    [bdepth, bbins] = meshgrid(depth, bins(1:end-1));
    bin_depth = bdepth + bbins;
    
    Z = bin_depth;
    ZmM = max(depth, [], 'omitnan');
    
    % Calculate shear per ensemble
    ensemble_shear_U = calc_ensemble_shear(U, bins);  % Implement this function
    ensemble_shear_V = calc_ensemble_shear(V, bins);  % Implement this function
    ensemble_shear_W = calc_ensemble_shear(W, bins);  % Implement this function
    
    % Create velocity dataframes for shear
    flatu = ensemble_shear_U(:);
    flatv = ensemble_shear_V(:);
    flatw = ensemble_shear_W(:);
    flatz = -Z(:);
    
    % Sort by depth
    [flatz, sort_ind] = sort(flatz, 'descend');
    flatu = flatu(sort_ind);
    flatv = flatv(sort_ind);
    flatw = flatw(sort_ind);
    
    % Shear binning
    [shear_binned, shear_binned_std, shear_cell_center, vels_in_bin] = bin_attr(flatu, flatv, flatw, flatz, dz, min(flatz, [], 'omitnan'));  % Implement this function
    
    % Shear to absolute
    [vel, vel_referenced, bin_centers] = shear_to_vel(shear_binned, shear_cell_center, [vx, vy, 0]);  % Implement this function
    
    % Define uncertainty of a single ping
    std_ping = 0.03;  % [m/s] (in average mode)
    vel_referenced_std = sqrt(std_ping^2 + shear_binned_std.^2);
end