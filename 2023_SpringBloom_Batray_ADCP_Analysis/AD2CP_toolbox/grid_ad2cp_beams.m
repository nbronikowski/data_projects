function ds = grid_ad2cp_beams(ds,Vrange)
    
    % Preallocate interpolated velocity matrices
    InterpVelocityBeam1 = nan(length(Vrange), length(ds.time));
    InterpVelocityBeam2 = nan(length(Vrange), length(ds.time));
    InterpVelocityBeam3 = nan(length(Vrange), length(ds.time));
    InterpVelocityBeam4 = nan(length(Vrange), length(ds.time));
    
    % Loop through each time (ping) and interpolate beam velocity onto the regular grid
    for x = 1:length(ds.time)
        InterpVelocityBeam1(:, x) = interp1(ds.TrueDepthBeam1(:, x), ds.VelocityBeam1(:, x), Vrange, 'linear', NaN);
        InterpVelocityBeam2(:, x) = interp1(ds.TrueDepthBeam2(:, x), ds.VelocityBeam2(:, x), Vrange, 'linear', NaN);
        InterpVelocityBeam3(:, x) = interp1(ds.TrueDepthBeam3(:, x), ds.VelocityBeam3(:, x), Vrange, 'linear', NaN);
        InterpVelocityBeam4(:, x) = interp1(ds.TrueDepthBeam4(:, x), ds.VelocityBeam4(:, x), Vrange, 'linear', NaN);
    end
    
    % Add the interpolated data back into the structure
    ds.InterpVelocityBeam1 = InterpVelocityBeam1;
    ds.InterpVelocityBeam2 = InterpVelocityBeam2;
    ds.InterpVelocityBeam3 = InterpVelocityBeam3;
    ds.InterpVelocityBeam4 = InterpVelocityBeam4;
    [~,ds.TrueDepth] = meshgrid(ds.time,Vrange);
    
end