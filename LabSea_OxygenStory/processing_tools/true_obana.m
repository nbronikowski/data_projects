function [y, y_m, y_err] = true_obana(x, posx, posy, gridx, gridy, r1x, r2x, r1y, r2y, method, signal_noise_ratio, n, nr)
    if nargin<13, nr=100; end
    if nargin<12, n=500; end
    if nargin<11, signal_noise_ratio=2; end
    if nargin<10, method=2; end
    if nargin<9, r2y=r2x; end
    if nargin<8, r1y=r1x; end
    
    sa = numel(gridx(:));
    signal_noise_ratio = 1/signal_noise_ratio + 1;
    
    % Initialize outputs
    y = NaN(sa, 1);
    y_err = NaN(sa, 1);
    y_m = NaN(sa, 1);
    
    % Pre-compute distances and indices for each grid point
    distances = cell(sa, 1); % Store indices of data points within influence radius
    for i = 1:sa
        Cx = abs(posx - gridx(i));
        Cy = abs(posy - gridy(i));
        distances{i} = find(Cx <= r2x & Cy <= r2y);
    end
    
    % Use parfor for parallel processing
    parfor istep = 1:sa
        gef = distances{istep};
        
        if isempty(gef)
            continue;
        elseif numel(gef) == 1
            y(istep) = x(gef);
            continue;
        end
        
        if numel(gef) > n
            gef = gef(randperm(numel(gef), n));
        end
        
        % Local covariance matrix calculation
        Ex_local = (posx(gef) - posx(gef)') / r1x;
        Ey_local = (posy(gef) - posy(gef)') / r1y;
        
        switch method
            case 2
                E_local = exp(-(Ex_local.^2 + Ey_local.^2));
            case 1
                E_local = exp(-sqrt(Ex_local.^2 + Ey_local.^2));
            otherwise
                error('Error: no valid model function');
        end
        
        E_local = E_local + diag(ones(size(gef)) / signal_noise_ratio);
        
        % Solve for weights and compute interpolated value
        weights = E_local \ ones(size(gef));
        weights = weights / sum(weights);
        
        y_m(istep) = sum(x(gef) .* weights);
        y(istep) = y_m(istep);
        
        % Simplified error estimation
        y_err(istep) = sqrt(1 - sum(weights .* (E_local * weights)));
    end
    
    y = reshape(y, size(gridx));
    y_err = reshape(y_err, size(gridx));
    y_m = reshape(y_m, size(gridx));
end



% function [y y_m y_err]=true_obana(x,posx,posy,gridx,gridy,r1x,r2x,r1y,r2y, method, signal_noise_ratio,n, nr)
%     if nargin<13, nr=100; % use 50 rand if not specified % use gauss if not specified
%         if nargin<12, n=500; % use 350 profile if not specified
%             if nargin<11, signal_noise_ratio=2;
%                 if nargin<10, method=2;
%                     if nargin<9, r2y=r2x; % make isotropic if not specified
%                         if nargin<8, r1y=r1x;
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     sn = numel(x);
%     sa = numel(gridx(:));
%     signal_noise_ratio = 1/signal_noise_ratio + 1;
%     
%     % Initialize outputs
%     y = NaN(sa, 1);
%     y_err = NaN(sa, 1);
%     y_m = NaN(sa, 1);
%     
%     % Loop over each grid point
%     for istep = 1:sa
%         % Find data points within the influence radius for the current grid point
%         Cx = abs(posx - gridx(istep));
%         Cy = abs(posy - gridy(istep));
%         gef = find(Cx <= r2x & Cy <= r2y);
%         
%         if isempty(gef)
%             continue; % No data points within the influence radius
%         elseif numel(gef) == 1
%             y(istep) = x(gef); % Only one data point within the influence radius
%             continue;
%         end
%         
%         % Select a random subset if the number of points is too large
%         if numel(gef) > n
%             gef = gef(randperm(numel(gef), n));
%         end
%         
%         % Local covariance matrix calculation
%         Ex_local = (posx(gef) - posx(gef)') / r1x;
%         Ey_local = (posy(gef) - posy(gef)') / r1y;
%         
%         switch method
%             case 2
%                 E_local = exp(-(Ex_local.^2 + Ey_local.^2));
%             case 1
%                 E_local = exp(-sqrt(Ex_local.^2 + Ey_local.^2));
%             otherwise
%                 error('Error: no valid model function');
%         end
%         
%         % Add the signal-to-noise ratio to the diagonal
%         E_local = E_local + diag(ones(size(gef)) / signal_noise_ratio);
% 
%         % Solve the system for weights
%         weights = E_local \ ones(size(gef));
%         weights = weights / sum(weights); % Normalize weights
%         
%         % Compute the interpolated value as the weighted mean
%         y_m(istep) = sum(x(gef) .* weights);
%         y(istep) = y_m(istep);
%         
%         % Error estimation
%         % In kriging, the estimation variance (kriging variance) is often used as a measure of uncertainty
%         % The kriging variance is the variance of the estimation error and can be computed as follows
%         % However, due to the simplification of using a local covariance matrix, the error estimation
%         % may not fully represent the true kriging variance as it would in a global kriging approach
%         y_err(istep) = sqrt(1 - sum(weights .* E_local * weights));
%     end
% 
%     % Reshape the output to the grid size
%     y = reshape(y, size(gridx));
%     y_err = reshape(y_err, size(gridx));
%     y_m = reshape(y_m, size(gridx));
% end