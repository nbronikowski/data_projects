function ds = beam2enu(ds,beam2xyz)
    % BEAM2ENU Convert velocity data from beam coordinates to ENU.
    % 
    % This function transforms velocity data from beam coordinates to XYZ to ENU.
    % Detailed explanation goes here
    
    % Transformation matrix layout
    % Beam: 1.   2.   3.   4.
    % X     1X.  2X.  3X.  4X.
    % Y     1Y.  2Y.  3Y.  4Y.
    % Z1    1Z1. 2Z1. 3Z1. 4Z1.
    % Z2    1Z2. 2Z2. 3Z2. 4Z2.
%     beam2xyz = ds.beam2xyz;

    % Assuming ds.attrs.burst_beam2xyz is a 4x4 matrix
    ds.AHRSRotationMatrix=calcAHRS(ds.heading,ds.roll,ds.pitch);

    % Initialize UVelocity, VVelocity, and WVelocity fields in ds
    ds.UVelocity = nan(size(ds.InterpVelocityBeam1));
    ds.VVelocity = nan(size(ds.InterpVelocityBeam1));
    ds.WVelocity = nan(size(ds.InterpVelocityBeam1));
    
    % Loop through each time point
    for x = 1:length(ds.time)
        if ds.pitch(x) < 0
            tot_vel = [ds.InterpVelocityBeam1(:,x), ds.InterpVelocityBeam2(:,x), ds.InterpVelocityBeam4(:,x)];
            beam2xyz_mat = beam2xyz(1:3, [1,2,4]);
        elseif ds.pitch(x) > 0
            tot_vel = [ds.InterpVelocityBeam2(:,x), ds.InterpVelocityBeam3(:,x), ds.InterpVelocityBeam4(:,x)];
            beam2xyz_mat = beam2xyz(1:3, 2:4);
        else
            tot_vel = [ds.VelocityBeam2(:,x), ds.VelocityBeam3(:,x), ds.VelocityBeam4(:,x)];
            beam2xyz_mat = beam2xyz(1:3, 2:4);
        end
        
        % If instrument is pointing down, bit 0 in status is equal to 1, rows 2 and 3 must change sign.
        beam2xyz_mat(2,:) = -beam2xyz_mat(2,:);
        beam2xyz_mat(3,:) = -beam2xyz_mat(3,:);
        
        % Convert to XYZ
        xyz = beam2xyz_mat* tot_vel';
        
        % Grab AHRS rotation matrix for this ping
        xyz2enuAHRS = reshape(ds.AHRSRotationMatrix(:,x), [3,3]);
        
        % Convert XYZ velocities to ENU
        enu = xyz2enuAHRS* xyz;
        
        % Assign to output structure
        ds.UVelocity(:,x) = enu(1,:)';
        ds.VVelocity(:,x) = enu(2,:)';
        ds.WVelocity(:,x) = enu(3,:)';
    end
end