function [O_ls, G_ls, bin_new, obs_per_bin] = ad2cp_inversion(U, V, dz, u_daverage, v_daverage, bins, depth, wDAC, wSmoothness)
    % INVERSION Perform inversion to separate ocean and glider velocities
    % from ADCP measurements.
    %
    % Inputs:
    % U - measured east-west velocities from ADCP
    % V - measured north-south velocities from ADCP
    % dz - desired vertical resolution
    % u_daverage - depth averaged velocity in U
    % v_daverage - depth averaged velocity in V
    % bins - bin depths for U and V measurements
    % depth - depth of the glider measured by the ADCP
    % wDAC - weight of the DAC constraint
    % wSmoothness - weight of the curvature minimizing constraint
    %
    % Outputs:
    % O_ls - ocean velocity profile
    % G_ls - glider velocity profile
    % bin_new - bin centers for the point in the profiles
    % obs_per_bin - number of good velocity observations per final profile bin

    %% Remove NaN rows and columns
    nanind = find(all(isnan(U), 2));
    if ~isempty(nanind)
        U(nanind, :) = [];
        V(nanind, :) = [];
        bins(nanind) = [];
    end
    
    nanind = find(all(isnan(U), 1));
    if ~isempty(nanind)
        U(:, nanind) = [];
        V(:, nanind) = [];
        depth(nanind) = [];
    end
    
    % Compute bin size and number
    bin_size = diff(bins(1:2));  % Assuming uniform bin size
    bin_num = length(bins);
    
    %% Compute actual depths of the ADCP bins
    [bdepth, bbins] = meshgrid(depth, bins);
    bin_depth = bdepth + bbins;
    Z = bin_depth;
    ZmM = max(depth, [], 'omitnan');
    
    %% Knowns from Visbeck (2002)
    nbin = size(U, 1);
    nt = size(U, 2);
    nd = nbin * nt;
    
    % Define bin edges and check each bin for data
    bin_edges = 0:dz:floor(max(bin_depth, [], 'all'));
    bin_count = NaN(length(bin_edges)-1, 1);
    
    for k = 1:(length(bin_edges)-1)
        ii = (bin_depth > bin_edges(k)) & (bin_depth <= bin_edges(k+1));
        bin_count(k) = sum(ii, 'all');
    end
    
    % Create list of bin centers
    bin_new = bin_edges(1:end-1) + dz/2;
    
%     % Adjust bins based on data availability
%     depth_ind = sum(bin_new < ZmM);
%     ind = find(bin_count > 0, 1, 'first');
%     bin_new = bin_new(ind:end);
%     z1 = bin_new(1);

    depth_ind = sum(bin_new > ZmM); % Adjusted to be equivalent to Python logic
    % Chop off the top of profile if no data
    ind = find(bin_count > 0, 1, 'first'); % Equivalent to np.argmax(bin_count > 0)
    bin_new = bin_new(ind:end); % Equivalent to Python slicing
    z1 = bin_new(1); % Equivalent to Python logic

    %% Create and populate G
    nz = length(bin_new);  % number of ocean velocities desired in output profile
    nm = nt + nz;          % G dimension (2), number of unknowns
    
    % Let's build the corresponding coefficient matrix G
    G = sparse(nd, nm);
    
    % Indexing of the G matrix was taken from Todd et al. 2012
    for ii = 1:nt           % Number of ADCP ensembles per segment
        for jj = 1:nbin     % Number of measured bins per ensemble
            
            % Uctd part of matrix
            G(nbin*(ii-1)+jj, ii) = -1;
            
            % This will fill in the Uocean part of the matrix. It loops through
            % all Z members and places them in the proper location in the G matrix
            % Find the difference between all bin centers and the current Z value
            dx = abs(bin_new - Z(jj, ii));
            
            % Find the minimum of these differences
            minx = nanmin(dx);
            
            % Finds bin_new index of the first match of Z and bin_new
            [~, idx] = min(dx-minx);
            
            % Uocean part of matrix
            G(nbin*(ii-1)+jj, nt+idx) = 1;
        end
    end
    
    %% Reshape U and V into the format of the d column vector
    % Based on how G is made, d needs to be ensembles stacked on one another vertically
    d_u = U(:);
    d_v = V(:);
    
    %% This chunk of code contains the constraints for depth averaged currents
    % Make sure the constraint is only applied to the final ocean velocity bins that the glider dives through
    % Don't apply it to the first bin and don't apply it to the bins below the gliders dive depth
    constraint = [zeros(1, nt), 0, repmat(dz, 1, nz-(1+depth_ind)), zeros(1, depth_ind)];
    
    % Ensure the L^2 norm of the constraint equation is unity
    constraint_norm = norm(constraint/ZmM);
    C = 1/constraint_norm;
    constraint_normalized = (C/ZmM) * constraint; % This is now equal to 1 (unity)
    
    % Build Gstar and add weight from Todd et al. 2017
    % Some smarts would be to calculate signal to noise ratio first
    Gstar = [G; wDAC * constraint_normalized];
    
    % Add the constraint for the depth averaged velocity from Todd et al. (2017)
    du = [d_u; wDAC * C * u_daverage];
    dv = [d_v; wDAC * C * v_daverage];
    d = complex(du, dv);

    %% THIS removes all NaN elements of d AND Gstar so the inversion doesn't blow up with NaNs
    ind2 = find(isnan(d));
    d(ind2) = [];
    
    Gstar = delete_rows_csr(Gstar, ind2);
    
    %% Test adding depth for tracking bin location
    % d is ensembles stacked on one another vertically so same for Z
    Z_filt = Z(:);
    Z_filt(ind2) = [];
    Z_filt = [Z_filt; 0];
    
    %% Calculation the number of observations per bin
    obs_per_bin = nan(1, length(bin_new));
    
    for x = 1:nz
        rows_where_nt_not_equal_zero = find(Gstar(1:length(Z_filt), nt+x) > 0);
        obs_per_bin(x) = length(rows_where_nt_not_equal_zero);
    end
    
    %% If there is no data in the last bin, drop that from the G matrix, bin_new, and obs_per_bin
    if obs_per_bin(end) == 0
        Gstar = Gstar(:, 1:end-1);
        bin_new = bin_new(1:end-1);
        obs_per_bin = obs_per_bin(1:end-1);
        % Update nz and nt
        nz = length(bin_new);
        nt = size(Gstar, 2) - nz;
    end
    
    %% Smoothness constraint
    % Only do this if the smoothness constraint is set
    if wSmoothness > 0
        % Add a vector of zeros, the length of nz, twice to the bottom of the data column vector
        d = [d; zeros(nz, 1); zeros(nz, 1)];
        
        % Constraint on smoothing Uocean side of matrix
        smoothing_matrix_Uocean = spdiags(repmat([-1, 2, -1], nz, 1), 0:2, nz, nz);
        smoothing_matrix1 = [zeros(nz, nt), smoothing_matrix_Uocean];
        
        % Constraint on smoothing Uglider side of matrix
        smoothing_matrix_Uglider = spdiags(repmat([-1, 2, -1], nz, 1), 0:2, nz, nt);
        smoothing_matrix2 = [smoothing_matrix_Uglider, zeros(nz, nz)];
        
        Gstar = [Gstar; wSmoothness * smoothing_matrix1; wSmoothness * smoothing_matrix2];
    end

    %% Run the Least-Squares Inversion!    
    x = lsqr(Gstar, d,1e-6,500);
    
    O_ls = x(nt+1:end);
    G_ls = x(1:nt);

end

% Define a function to delete rows from sparse matrix
function mat = delete_rows_csr(mat, indices)
    % Remove the rows denoted by `indices` from the CSR sparse matrix `mat`.
    if ~isa(mat, 'double')
        error('works only for CSR format -- use .tocsr() first');
    end
    mask = true(1, size(mat, 1));
    mask(indices) = false;
    mat = mat(mask, :);
end
    