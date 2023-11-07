function corrected_heading=correct_AD2CP_heading(time,head,pitch,roll,pressure,x,y,z,ptcheck)
    headt=[];
    pitcht=[]; 
    rollt=[]; 
    pressuret=[]; 
    xt=[]; 
    yt=[]; 
    zt=[]; 
    timet=[];

    [head_M,head_N]=size(head);
  
    % Store in temporary variables
    headt=[headt head']; 
    pitcht=[pitcht pitch']; 
    rollt=[rollt roll']; 
    pressuret=[pressuret pressure']; 
    xt=[xt x']; 
    yt=[yt y']; 
    zt=[zt z']; 
    timet=[timet time'];
       
    % Rename and transpose the temporary variables
    head=headt'; 
    pitch=pitcht'; 
    roll=rollt'; 
    pressure=pressuret'; 
    x=xt'; 
    y=yt'; 
    z=zt';
    time=timet';
    
    clear headt pitcht rollt pressuret xt yt zt timet timeBT

    XYZ_original=[x y z];    
    %% Pitch Dependent magnetic heading correction
    % The pitch of the glider changes by moving the battery pack. The battery
    % pack produces a magnetic field so when it moves, the field moves and this
    % effects our compass data. Having all of the data from a mission loaded in
    % for this correction ensures the best fit possible.
    
    pitch_ranges=(-20:1:20);
    for k=1:length(pitch_ranges)-1
        % Find where pitch falls into this particular range   
        ii=find(pitch > pitch_ranges(k) & pitch < pitch_ranges(k+1));
        if ~isempty(ii)
            % Index the magnetic field variables to only be referencing the current
            % pitch range
            XYZ1=[x(ii) y(ii) z(ii)];
         
            % Fits the given cartesian coordinates to an ellipsoid and gives the
            % center of that shape as the first output.
            [offset, ~, ~, ~, ~ ] = ellipsoid_fit_tnm( XYZ1 );
         
            % Creates xyz variables that are offset by the center of the
            % magnetic field recorded in the data
            x1 = x(ii)-offset(1);
            y1 = y(ii)-offset(2);
            z1 = z(ii)-offset(3);
        
            % Find the center of the new magnetic data that should be closer 0,0,z
            % than the original data.
            [new_center, ~, ~, ~, ~ ] = ellipsoid_fit_tnm([x1 y1 z1]);
            
            % If the center of the offset magnetic data is not near 0,0 in the x
            % and y planes, we keep the original values.
            if abs(new_center(1)) > 150 || abs(new_center(2)) > 150
                x1=x(ii); y1=y(ii);z1=z(ii);
            end
        
            % Replace the indexed values with the new corrected values
            % If the previous if statement is true, nothing happens here, x(ii)=x1
            % which = x(ii).
            % If it is not, then x(ii)=x(ii)-offset1 etc.
            x(ii)=x1; 
            y(ii)=y1; 
            z(ii)=z1;
        end
        clear x1 y1 z1 offset radii evecs evals pars
    end
    
    XYZ_final=[x y z];
    
    %% Re-calculate heading with magnetometer corrections
    headingf = NaN*pressure; headingo=headingf;
    for k=1:length(pressure)
        % Final heading
        headingf(k) = CalcMidlifeHeading(XYZ_final(k,:), pitch(k), roll(k), 1);
        % Original heading
        headingo(k) = CalcMidlifeHeading(XYZ_original(k,:), pitch(k), roll(k), 1);
    end
    
    % Plots the corrected (blue) and uncorrected (red) heading values to check
    % the correction was successful.
    corrected_heading = reshape(headingf,[head_M,head_N]);

    if ptcheck>0
        figure(); plot(XYZ_original(:,1),XYZ_original(:,2),'.r',x,y,'.b')
        legend('uncorrected (red)','corrected (blue)')
        xlabel('magnetic x')
        ylabel('magnetic y')
        close gcf
    end
end