function RotMatrix = calcAHRS(headingVal,rollVal,pitchVal)
    
RotMatrix = NaN(9,length(pitchVal));
    % depending on pitch use different beams from matrix.
    for k = 1:length(pitchVal)

        % If instrument pointed down we need to flip sign on rows...        
        hh = deg2rad(headingVal(k)-90);
        pp = deg2rad(pitchVal(k));
        rr = deg2rad(rollVal(k));
        
        % Make heading matrix
        H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
        
        % Make tilt matrix
        P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
              0             cos(rr)          -sin(rr);  ...
              sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
        
        % Make resulting transformation matrix
        R = H*P;
        RotMatrix(:,k) = R(:);
    end
end